!> The main IO server functionality which handles waiting for commands and data both of which are delt with.
!! The lower level details of the communication, configuration parsing etc are all held elsewhere. The server
!! can be thought of similar to a bus, with command and data channels. The command gives context to what is on
!! the data channel and not all commands require data (such as deregistration of MONC process)
module io_server_mod
  use datadefn_mod, only : DEFAULT_PRECISION, STRING_LENGTH, LONG_STRING_LENGTH
  use configuration_parser_mod, only : DATA_SIZE_STRIDE, io_configuration_type, io_configuration_data_definition_type, &
       io_configuration_registered_monc_type, configuration_parse, extend_registered_moncs_array, retrieve_data_definition, &
       build_definition_description_type_from_configuration, build_field_description_type_from_configuration, get_monc_location, &
       get_io_xml
  use mpi_communication_mod, only : free_mpi_type, get_number_io_servers, get_my_io_rank, lock_mpi, unlock_mpi, &
       waitall_for_mpi_requests, initialise_mpi_communication
  use diagnostic_federator_mod, only : initialise_diagnostic_federator, finalise_diagnostic_federator, &
       pass_fields_to_diagnostics_federator, determine_diagnostics_fields_available, &
       diagnostics_federator_task
  use writer_federator_mod, only : initialise_writer_federator, finalise_writer_federator, check_writer_for_trigger, &
       inform_writer_federator_fields_present, inform_writer_federator_time_point, provide_q_field_names_to_writer_federator
  use writer_field_manager_mod, only : initialise_writer_field_manager, finalise_writer_field_manager, &
       provide_monc_data_to_writer_federator, writer_federator_task_from_monc
  use collections_mod, only : hashset_type, hashmap_type, map_type, iterator_type, c_get_integer, c_put_integer, c_is_empty, &
       c_remove, c_add_string, c_integer_at, c_free, c_get_iterator, c_has_next, c_next_mapentry
  use conversions_mod, only : conv_to_string
  use string_utils_mod, only : replace_character
  use optionsdatabase_mod, only : options_get_string
  use io_server_client_mod, only : REGISTER_COMMAND, DEREGISTER_COMMAND, INTER_IO_COMMUNICATION, DATA_COMMAND_START, DATA_TAG, &
       LOCAL_SIZES_KEY, LOCAL_START_POINTS_KEY, LOCAL_END_POINTS_KEY, NUMBER_Q_INDICIES_KEY, SCALAR_FIELD_TYPE, &
       data_sizing_description_type, definition_description_type, field_description_type,&
       get_data_description_from_name, build_mpi_type_field_description, build_mpi_type_definition_description
  use forthread_mod, only : forthread_init, forthread_rwlock_rdlock, forthread_rwlock_wrlock, forthread_rwlock_tryrdlock, &
       forthread_rwlock_unlock, forthread_rwlock_init, forthread_rwlock_destroy, forthread_mutex_init, forthread_mutex_lock, &
       forthread_mutex_unlock, forthread_cond_wait, forthread_cond_signal, forthread_cond_init
  use data_utils_mod, only : check_thread_status, build_depacking_description
  use logging_mod, only : LOG_ERROR, LOG_WARN, log_log, initialise_logging
  use mpi, only : MPI_COMM_WORLD, MPI_INT
  use io_server_state_reader_mod, only : read_io_server_configuration

  use edat
  use iso_c_binding, only : c_ptr, c_int
  implicit none

#ifndef TEST_MODE
  private
#endif

  integer :: mpi_type_data_sizing_description, & !< The MPI type for field sizing (i.e. array size etc send when MONCs register)
       mpi_type_definition_description, & !< The MPI data type for data descriptions sent to MONCs
       mpi_type_field_description !< The MPI data type for field descriptions sent to MONCs
  type(io_configuration_type), volatile, target, save :: io_configuration !< Internal representation of the IO configuration
  logical, volatile :: initialised_present_data
  type(field_description_type), dimension(:), allocatable :: registree_field_descriptions
  type(definition_description_type), dimension(:), allocatable :: registree_definition_descriptions

  integer, volatile :: monc_registration_lock

  public io_server_run
contains

  !> Called to start the IO server and once this subroutine returns then it indicates that the IO server has finished.
  !! The runtime is spent in here awaiting commands and then dealing with them. Termination occurs when all MONC processes
  !! have deregistered, note that to trigger this then at least one MONC process must first register
  !! @param io_communicator_arg The IO communicator containing just the IO servers
  !! @param io_xml_configuration Textual XML configuration that is used to set up the IO server
  subroutine io_server_run(options_database, io_communicator_arg, &
       provided_threading, total_global_processes, continuation_run, io_configuration_file)
    type(hashmap_type), intent(inout) :: options_database
    integer, intent(in) :: io_communicator_arg, provided_threading, total_global_processes
    logical, intent(in) :: continuation_run
    character(len=LONG_STRING_LENGTH), intent(in) :: io_configuration_file

    integer :: command, source, my_rank, ierr
    character, dimension(:), allocatable :: data_buffer, io_xml_configuration
    type(hashmap_type) :: diagnostic_generation_frequency

    if (continuation_run) then
      ! Handle case where we need to allocate this due to no IO server config
      call read_io_server_configuration(options_get_string(options_database, "checkpoint"), &
           io_xml_configuration, io_communicator_arg)
    end if

    if (.not. allocated(io_xml_configuration)) then
      io_xml_configuration=get_io_xml(io_configuration_file)
      if (continuation_run) then
        call mpi_comm_rank(io_communicator_arg, my_rank, ierr)
        if (my_rank == 0) then
          call log_log(LOG_WARN, "No IO server configuration in checkpoint file - starting from XML provided file instead")
        end if
      end if
    end if

    call configuration_parse(options_database, io_xml_configuration, io_configuration)
    deallocate(io_xml_configuration)

    call edatInitialiseWithCommunicator(io_communicator_arg)

    call check_thread_status(forthread_init())

    call initialise_mpi_communication(provided_threading)
    call check_thread_status(forthread_rwlock_init(monc_registration_lock, -1))
    call check_thread_status(forthread_mutex_init(io_configuration%general_info_mutex, -1))
    initialised_present_data=.false.
    io_configuration%io_communicator=io_communicator_arg
    io_configuration%number_of_io_servers=get_number_io_servers(io_communicator_arg)
    io_configuration%number_of_global_moncs=total_global_processes-io_configuration%number_of_io_servers
    io_configuration%my_io_rank=get_my_io_rank(io_communicator_arg)
    call initialise_logging(io_configuration%my_io_rank)
    registree_definition_descriptions=build_definition_description_type_from_configuration(io_configuration)
    registree_field_descriptions=build_field_description_type_from_configuration(io_configuration)
    diagnostic_generation_frequency=initialise_diagnostic_federator(io_configuration)
    call initialise_writer_federator(io_configuration, diagnostic_generation_frequency, continuation_run)
    call c_free(diagnostic_generation_frequency)
    call initialise_writer_field_manager(io_configuration, continuation_run)

    mpi_type_definition_description=build_mpi_type_definition_description()
    mpi_type_field_description=build_mpi_type_field_description()

    call edatSubmitPersistentTask(register_monc, 1, EDAT_ANY, "register")
    call edatSubmitNamedTask(register_monc, "placeholder", 1, 0, "placeholder")

    call edatFinalise()
    call free_individual_registered_monc_aspects()
    call finalise_writer_field_manager()
    call finalise_writer_federator()
    call finalise_diagnostic_federator(io_configuration)
    call check_thread_status(forthread_rwlock_destroy(monc_registration_lock))
    call free_mpi_type(mpi_type_definition_description)
    call free_mpi_type(mpi_type_field_description)
  end subroutine io_server_run

  !> Frees up the memory associated with individual registered MONCs. This is done at the end for all MONCs as we can't
  !! deallocate dynamically in a threaded environment without excessive ordering and locking in case some data processing
  !! is queued or in progress
  subroutine free_individual_registered_monc_aspects()
    integer :: i

    do i=1, size(io_configuration%registered_moncs)
      if (allocated(io_configuration%registered_moncs(i)%field_start_locations)) &
           deallocate(io_configuration%registered_moncs(i)%field_start_locations)
      if (allocated(io_configuration%registered_moncs(i)%field_end_locations)) &
           deallocate(io_configuration%registered_moncs(i)%field_end_locations)
      if (allocated(io_configuration%registered_moncs(i)%definition_names)) &
           deallocate(io_configuration%registered_moncs(i)%definition_names)
      if (allocated(io_configuration%registered_moncs(i)%dimensions)) deallocate(io_configuration%registered_moncs(i)%dimensions)
    end do
  end subroutine free_individual_registered_monc_aspects

  subroutine handle_data(events, number_events)
	  type(c_ptr), intent(in), target :: events
    integer(c_int), value, intent(in) :: number_events

    type(EDAT_Event) :: processed_events(number_events)
    integer :: monc_location, data_set, source, matched_datadefn_index, ierr

    character(len=100) :: base_task_eid

	  call getEvents(events, number_events,  processed_events)

	  source=processed_events(1)%metadata%source
    data_set=processed_events(1)%int_data(1)

    call check_thread_status(forthread_rwlock_rdlock(monc_registration_lock))
    monc_location=get_monc_location(io_configuration, source)

    matched_datadefn_index=retrieve_data_definition(io_configuration, &
         io_configuration%registered_moncs(monc_location)%definition_names(data_set))

    call check_thread_status(forthread_rwlock_unlock(monc_registration_lock))

    if (matched_datadefn_index .gt. 0) then
      call inform_writer_federator_time_point(io_configuration, source, data_set, processed_events(2)%byte_data)
      base_task_eid=trim(conv_to_string(source)) // "#" // trim(conv_to_string(data_set))
      call edatSubmitTask(diagnostics_federator_task, 1, EDAT_SELF, trim(base_task_eid)//"_diag")
      call edatFireEvent(processed_events(2)%byte_data, EDAT_BYTE, processed_events(2)%metadata%number_elements, &
        EDAT_SELF, trim(base_task_eid)//"_diag")

      call edatSubmitTask(writer_federator_task_from_monc, 1, EDAT_SELF, trim(base_task_eid)//"_writer")
      call edatFireEvent(processed_events(2)%byte_data, EDAT_BYTE, processed_events(2)%metadata%number_elements, &
        EDAT_SELF, trim(base_task_eid)//"_writer")

      !call pass_fields_to_diagnostics_federator(io_configuration, source, data_set, processed_events(2)%byte_data)
      !call provide_monc_data_to_writer_federator(io_configuration, source, data_set, processed_events(2)%byte_data)
      call check_writer_for_trigger(io_configuration, source, data_set, processed_events(2)%byte_data)
    else
      call log_log(LOG_WARN, "IO server can not find matching data definition with name "&
           //io_configuration%registered_moncs(monc_location)%definition_names(data_set))
    end if
		!call lock_mpi()
		!call MPI_Send(0, 0, MPI_INT, source, 0, MPI_COMM_WORLD, ierr)
		!call unlock_mpi()
  end subroutine handle_data

  subroutine register_monc(events, number_events)
	  type(c_ptr), intent(in), target :: events
    integer(c_int), value, intent(in) :: number_events

    integer :: configuration_send_request(2), ierr, number_data_definitions, this_monc_index, source

    type(EDAT_Event) :: processed_events(number_events)
	  call getEvents(events, number_events,  processed_events)

    source=processed_events(1)%metadata%source
    configuration_send_request=send_configuration_to_registree(source)
    number_data_definitions=io_configuration%number_of_data_definitions

    call check_thread_status(forthread_rwlock_wrlock(monc_registration_lock))

    io_configuration%number_of_moncs=io_configuration%number_of_moncs+1
    this_monc_index=io_configuration%number_of_moncs
    if (io_configuration%number_of_moncs .gt. size(io_configuration%registered_moncs)) then
      call log_log(LOG_ERROR, "You have a high ratio of computational cores to IO servers, the limit is currently 100")
      ! The extension of the MONC registration array is broken as the pointers involved in the map does not get copied across
      ! we could manually do this, but that is for another day! If you need to extend these limits either increase the constants
      ! or fix the extension, I don't think it will be too hard to fix the extension bit (copy the maps manually)
      call extend_registered_moncs_array(io_configuration)
    end if

    io_configuration%active_moncs=io_configuration%active_moncs+1
    call check_thread_status(forthread_rwlock_unlock(monc_registration_lock))

    call c_put_integer(io_configuration%monc_to_index, conv_to_string(source), this_monc_index)

    call check_thread_status(forthread_mutex_init(io_configuration%registered_moncs(this_monc_index)%active_mutex, -1))
    call check_thread_status(forthread_cond_init(&
         io_configuration%registered_moncs(this_monc_index)%deactivate_condition_variable, -1))
    io_configuration%registered_moncs(this_monc_index)%active_threads=0
    io_configuration%registered_moncs(this_monc_index)%source_id=source

    allocate(io_configuration%registered_moncs(this_monc_index)%field_start_locations(number_data_definitions), &
         io_configuration%registered_moncs(this_monc_index)%field_end_locations(number_data_definitions), &
         io_configuration%registered_moncs(this_monc_index)%definition_names(number_data_definitions), &
         io_configuration%registered_moncs(this_monc_index)%dimensions(number_data_definitions))

    ! Wait for configuration to have been sent to registree
    call waitall_for_mpi_requests(configuration_send_request, 2)

    call edatSubmitTask(monc_registration_handshake, 2, source, "fieldsizes", source, "moncinfo")
  end subroutine register_monc

  subroutine monc_registration_handshake(events, number_events)
	  type(c_ptr), intent(in), target :: events
    integer(c_int), value, intent(in) :: number_events

    type(EDAT_Event) :: processed_events(number_events)
    type(io_configuration_registered_monc_type), pointer :: monc_defn
    integer :: monc_location, depacking_description, i, source
    type(data_sizing_description_type) :: field_description
    logical :: field_found, desched
    type(data_sizing_description_type) :: q, data_description(io_configuration%number_of_distinct_data_fields+4)

	  call getEvents(events, number_events,  processed_events)

		source=processed_events(1)%metadata%source

    data_description=transfer(processed_events(1)%byte_data, data_description)
    monc_location=get_monc_location(io_configuration, source)
    monc_defn=>io_configuration%registered_moncs(monc_location)
    call handle_monc_dimension_information(data_description, monc_defn)

    do i=1, io_configuration%number_of_data_definitions
      call build_depacking_description(io_configuration%data_definitions(i), data_description, &
        monc_defn%field_start_locations(i), monc_defn%field_end_locations(i), monc_defn%dimensions(i))

      monc_defn%definition_names(i)=io_configuration%data_definitions(i)%name
    end do
    if (.not. initialised_present_data) then
      initialised_present_data=.true.
      field_found=get_data_description_from_name(data_description, NUMBER_Q_INDICIES_KEY, field_description)
      call c_put_integer(io_configuration%dimension_sizing, "active_q_indicies", field_description%dim_sizes(1))
      call register_present_field_names_to_federators(data_description, size(data_description))
    end if
    call get_monc_information_data(processed_events(2)%byte_data)
		call edatSubmitTask(deregister, 1, source, "deregister")
    call edatSubmitPersistentTask(handle_data, 2, source, "sendreq", source, "data")
    desched=edatRemoveTask("placeholder")
  end subroutine monc_registration_handshake

  subroutine deregister(events, number_events)
	  type(c_ptr), intent(in), target :: events
    integer(c_int), value, intent(in) :: number_events

    integer :: monc_location
    type(EDAT_Event) :: processed_events(number_events)

	  call getEvents(events, number_events,  processed_events)

		return

    monc_location=get_monc_location(io_configuration, processed_events(1)%metadata%source)
    call check_thread_status(forthread_mutex_lock(io_configuration%registered_moncs(monc_location)%active_mutex))
    do while (io_configuration%registered_moncs(monc_location)%active_threads .gt. 0)
      call check_thread_status(forthread_cond_wait(io_configuration%registered_moncs(monc_location)%deactivate_condition_variable,&
             io_configuration%registered_moncs(monc_location)%active_mutex))
    end do
    call check_thread_status(forthread_mutex_unlock(io_configuration%registered_moncs(monc_location)%active_mutex))
    call check_thread_status(forthread_rwlock_wrlock(monc_registration_lock))
    io_configuration%active_moncs=io_configuration%active_moncs-1
    call check_thread_status(forthread_rwlock_unlock(monc_registration_lock))
  end subroutine deregister

  !> Sends the data and field descriptions to the MONC process that just registered with the IO server
  !! @param source The MPI rank (MPI_COMM_WORLD) of the registree
  !! @returns The nonblocking send request handles which can be waited for completion later (overlap compute and communication)
  function send_configuration_to_registree(source)
    integer, intent(in) :: source
    integer :: send_configuration_to_registree(2)

    integer :: ierr, srequest(2)

    call lock_mpi()
    call mpi_isend(registree_definition_descriptions, size(registree_definition_descriptions), mpi_type_definition_description, &
         source, DATA_TAG, MPI_COMM_WORLD, srequest(1), ierr)
    call mpi_isend(registree_field_descriptions, size(registree_field_descriptions), mpi_type_field_description, &
         source, DATA_TAG, MPI_COMM_WORLD, srequest(2), ierr)
    call unlock_mpi()

    send_configuration_to_registree=srequest
  end function send_configuration_to_registree

  !> Retrieves MONC information data, this is sent by MONC (and received) regardless, but only actioned if the data has not
  !! already been set
  !! @param source MONC source process
  subroutine get_monc_information_data(buffer)
    character, dimension(:), intent(in) :: buffer

    character(len=STRING_LENGTH) :: q_field_name
    integer :: z_size, num_q_fields, n, current_point
    type(data_sizing_description_type) :: field_description
    real(kind=DEFAULT_PRECISION) :: dreal
    logical :: field_found

    z_size=c_get_integer(io_configuration%dimension_sizing, "z")
    num_q_fields=c_get_integer(io_configuration%dimension_sizing, "qfields")

    if (.not. io_configuration%general_info_set) then
      call check_thread_status(forthread_mutex_lock(io_configuration%general_info_mutex))
      if (.not. io_configuration%general_info_set) then
        io_configuration%general_info_set=.true.
        allocate(io_configuration%zn_field(z_size))
        io_configuration%zn_field=transfer(buffer(1:kind(dreal)*z_size), io_configuration%zn_field)
        current_point=(kind(dreal)*z_size)
        if (num_q_fields .gt. 0) then
          do n=1, num_q_fields
            q_field_name=transfer(buffer(current_point+1:current_point+STRING_LENGTH), q_field_name)
            current_point=current_point+STRING_LENGTH
            call replace_character(q_field_name, " ", "_")
            call c_add_string(io_configuration%q_field_names, q_field_name)
          end do
        end if
      end if
      call provide_q_field_names_to_writer_federator(io_configuration%q_field_names)
      call check_thread_status(forthread_mutex_unlock(io_configuration%general_info_mutex))
    end if
  end subroutine get_monc_information_data

  !> Registers with the writer federator the set of fields (prognostic and diagnostic) that are available, this is based on
  !! the array/optional fields present from MONC and the non-optional scalars. This is quite an expensive operation, so only
  !! done once
  !! @param data_description Array of data descriptions from MONC
  !! @param recv_count Number of data descriptions
  subroutine register_present_field_names_to_federators(data_description, recv_count)
    type(data_sizing_description_type), dimension(:), intent(in) :: data_description
    integer, intent(in) :: recv_count

    type(hashset_type) :: present_field_names
    type(hashmap_type) :: diagnostics_field_names_and_roots
    integer :: i, j

    do i=1, recv_count
      call c_add_string(present_field_names, data_description(i)%field_name)
    end do
    do i=1, io_configuration%number_of_data_definitions
      do j=1, io_configuration%data_definitions(i)%number_of_data_fields
        if (io_configuration%data_definitions(i)%fields(j)%field_type == SCALAR_FIELD_TYPE .and. .not. &
             io_configuration%data_definitions(i)%fields(j)%optional) then
          call c_add_string(present_field_names, io_configuration%data_definitions(i)%fields(j)%name)
        end if
      end do
    end do
    call c_add_string(present_field_names, "time")
    call c_add_string(present_field_names, "timestep")
    call inform_writer_federator_fields_present(io_configuration, present_field_names)
    diagnostics_field_names_and_roots=determine_diagnostics_fields_available(present_field_names)
    call inform_writer_federator_fields_present(io_configuration, diag_field_names_and_roots=diagnostics_field_names_and_roots)
    call c_free(present_field_names)
    call c_free(diagnostics_field_names_and_roots)
  end subroutine register_present_field_names_to_federators

  !> Handles the provided local MONC dimension and data layout information
  !! @param data_description The data descriptions sent over from MONC
  !! @param monc_defn The corresponding MONC definition data structure
  subroutine handle_monc_dimension_information(data_description, monc_defn)
    type(io_configuration_registered_monc_type), intent(inout) :: monc_defn
    type(data_sizing_description_type), dimension(:) :: data_description

    type(data_sizing_description_type) :: field_description
    integer :: i
    logical :: field_found

    field_found=get_data_description_from_name(data_description, LOCAL_SIZES_KEY, field_description)
    if (.not. field_found) call log_log(LOG_ERROR, "Malformed MONC registration, no local size information")
    do i=1,3
      monc_defn%local_dim_sizes(i)=field_description%dim_sizes(i)
    end do
    field_found=get_data_description_from_name(data_description, LOCAL_START_POINTS_KEY, field_description)
    if (.not. field_found) call log_log(LOG_ERROR, "Malformed MONC registration, no local start point information")
    do i=1,3
      monc_defn%local_dim_starts(i)=field_description%dim_sizes(i)
    end do
    field_found=get_data_description_from_name(data_description, LOCAL_END_POINTS_KEY, field_description)
    if (.not. field_found) call log_log(LOG_ERROR, "Malformed MONC registration, no local end point information")
    do i=1,3
      monc_defn%local_dim_ends(i)=field_description%dim_sizes(i)
    end do
  end subroutine handle_monc_dimension_information
end module io_server_mod
