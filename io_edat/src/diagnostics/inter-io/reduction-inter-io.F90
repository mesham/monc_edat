!> Reduction inter IO action which will perform reductions between IO servers. This is not as trivial as calling the MPI function
!! as it is nondeterministic when messages will arrive and hence when one reduction on a process and a reduction on another
!! should be called.
module reduction_inter_io_mod
  use datadefn_mod, only : DEFAULT_PRECISION, DOUBLE_PRECISION, STRING_LENGTH
  use configuration_parser_mod, only : io_configuration_type, io_configuration_inter_communication_description
  use collections_mod, only : map_type, hashmap_type, list_type, iterator_type, mapentry_type, c_get_generic, c_get_string, &
       c_add_string, c_put_generic, c_remove, c_is_empty, c_contains, c_free, c_generic_at, c_get_iterator, c_has_next, &
       c_next_mapentry, c_next_string
  use conversions_mod, only : conv_to_string, conv_to_integer
  use forthread_mod, only : forthread_mutex_init, forthread_mutex_lock, forthread_mutex_trylock, &
       forthread_mutex_unlock, forthread_mutex_destroy, forthread_rwlock_rdlock, forthread_rwlock_wrlock, &
       forthread_rwlock_unlock, forthread_rwlock_init, forthread_rwlock_destroy, forthread_rwlock_trywrlock
  use data_utils_mod, only : check_thread_status
  use logging_mod, only : LOG_ERROR, log_log
  use mpi, only : MPI_DOUBLE_PRECISION, MPI_INT, MPI_ANY_SOURCE, MPI_REQUEST_NULL, MPI_STATUS_IGNORE, MPI_CHARACTER, MPI_BYTE
  use mpi_communication_mod, only : lock_mpi, unlock_mpi, wait_for_mpi_request

  use edat
  use iso_c_binding, only : c_ptr, c_int
  implicit none

#ifndef TEST_MODE
  private
#endif

  abstract interface
     subroutine handle_completion(io_configuration, values, field_name, timestep)
       import DEFAULT_PRECISION, STRING_LENGTH, io_configuration_type
       type(io_configuration_type), intent(inout) :: io_configuration
       real(DEFAULT_PRECISION), dimension(:) :: values
       character(len=STRING_LENGTH) :: field_name
       integer :: timestep
     end subroutine handle_completion
  end interface

  !< The types of reduction operator that are supported
  integer, parameter :: MEAN=1, MIN=2, MAX=3, SUM=4

  type reduction_progress_type
     real(DEFAULT_PRECISION), dimension(:), allocatable :: values
     character(len=STRING_LENGTH) :: field_name
     integer :: contributed_moncs, timestep, reduction_operator, mutex, root
     procedure(handle_completion), pointer, nopass :: completion_procedure
  end type reduction_progress_type

  integer, volatile :: reduction_progress_rwlock
  type(hashmap_type), volatile :: reduction_progresses
  logical, volatile :: initialised=.false.

  type(io_configuration_type), pointer :: this_io_configuration

  public init_reduction_inter_io, finalise_reduction_inter_io, perform_inter_io_reduction, &
    perform_inter_io_allreduction, get_reduction_operator
contains

  !> Initialises the reduction action
  !! @param io_configuration The IO server configuration
  subroutine init_reduction_inter_io(io_configuration)
    type(io_configuration_type), target, intent(inout) :: io_configuration

    this_io_configuration=>io_configuration

    if (.not. initialised) then
      initialised=.true.
      call check_thread_status(forthread_rwlock_init(reduction_progress_rwlock, -1))
    end if
  end subroutine init_reduction_inter_io

  !> Finalises the reduction action, waiting for all outstanding requests and then freeing data
  !! @param io_configuration Configuration state of the IO server
  subroutine finalise_reduction_inter_io(io_configuration)
    type(io_configuration_type), intent(inout) :: io_configuration

    type(reduction_progress_type) :: progress
    type(iterator_type) :: iterator

    if (initialised) then
      call check_thread_status(forthread_rwlock_destroy(reduction_progress_rwlock))
      initialised=.false.
    end if
  end subroutine finalise_reduction_inter_io

  subroutine perform_inter_io_allreduction(io_configuration, field_values, field_size, reduction_field_name, reduction_op, &
      timestep, completion_procedure)
    type(io_configuration_type), intent(inout) :: io_configuration
    real(kind=DOUBLE_PRECISION), dimension(:) :: field_values
    integer, intent(in) :: field_size, reduction_op, timestep
    character(len=*), intent(in) :: reduction_field_name
    procedure(handle_completion) :: completion_procedure

    call do_inter_io_reduction(io_configuration, field_values, field_size, reduction_field_name, reduction_op, &
       .false., -1, timestep, completion_procedure)
  end subroutine perform_inter_io_allreduction

  subroutine perform_inter_io_reduction(io_configuration, field_values, field_size, reduction_field_name, reduction_op, &
       root, timestep, completion_procedure)
    type(io_configuration_type), intent(inout) :: io_configuration
    real(kind=DOUBLE_PRECISION), dimension(:) :: field_values
    integer, intent(in) :: field_size, reduction_op, root, timestep
    character(len=*), intent(in) :: reduction_field_name
    procedure(handle_completion) :: completion_procedure

    call do_inter_io_reduction(io_configuration, field_values, field_size, reduction_field_name, reduction_op, &
       .true., root, timestep, completion_procedure)
  end subroutine perform_inter_io_reduction

  !> Actually handles the processing for this data wrt the vertical reduction
  !! @param io_configuration Configuration of the IO server
  !! @param field_values The values to communicate
  !! @param field_size Number of elements to communicate
  !! @param reduction_field_name Field name that the reduction will be performed over
  !! @param reduction_op The reduction operator to use
  !! @param root The root IO server process
  !! @param timestep The timestep this is issued at
  !! @param completion_procedure Callback completion procedure
  subroutine do_inter_io_reduction(io_configuration, field_values, field_size, reduction_field_name, reduction_op, &
       has_root, root, timestep, completion_procedure)
    type(io_configuration_type), intent(inout) :: io_configuration
    real(kind=DOUBLE_PRECISION), dimension(:) :: field_values
    integer, intent(in) :: field_size, reduction_op, root, timestep
    logical, intent(in) :: has_root
    character(len=*), intent(in) :: reduction_field_name
    procedure(handle_completion) :: completion_procedure

    type(reduction_progress_type), pointer :: reduction_progress
    logical :: collective_values_new

    reduction_progress=>find_or_add_reduction_progress(timestep, reduction_op, has_root, root, io_configuration%my_io_rank, &
      reduction_field_name, completion_procedure)

    call check_thread_status(forthread_mutex_lock(reduction_progress%mutex))
    reduction_progress%contributed_moncs=reduction_progress%contributed_moncs+1

    collective_values_new=.not. allocated(reduction_progress%values)
    if (collective_values_new) allocate(reduction_progress%values(field_size))

    call integrate_io_server_collective_values(reduction_op, reduction_progress, field_values, field_size, collective_values_new)
    if (reduction_progress%contributed_moncs == io_configuration%number_of_moncs) then
      call check_thread_status(forthread_mutex_unlock(reduction_progress%mutex))
      call handle_local_moncs_completed_collective(io_configuration, has_root, reduction_progress)
    else
      call check_thread_status(forthread_mutex_unlock(reduction_progress%mutex))
    end if
  end subroutine do_inter_io_reduction

  !> Integrates the collective values from another IO server into the currently stored values
  !! @param reduction_op The reduction operator to perform
  !! @param reduction_progress The progress data type which is updated
  !! @param single_server_values The values from the IO server which need to be integrated
  !! @param dim_one_size Size in first dimension
  !! @param target_index The index where we are putting the values into the current value array of reduction progress
  !! @param collective_values_empty Whether the collective values is empty
  subroutine integrate_io_server_collective_values(reduction_op, reduction_progress, single_server_values, &
       number_elements, collective_values_empty)
    integer, intent(in) :: reduction_op, number_elements
    logical, intent(in) :: collective_values_empty
    type(reduction_progress_type), intent(inout) :: reduction_progress
    real(kind=DOUBLE_PRECISION), dimension(:), intent(in) :: single_server_values

    integer :: k

    if (collective_values_empty) then
      reduction_progress%values=single_server_values
    else
      if (reduction_op == MEAN .or. reduction_op == SUM) then
        reduction_progress%values=reduction_progress%values+single_server_values
      else if (reduction_op == MIN .or. reduction_op == MAX) then
        do k=1, number_elements
          if (reduction_op == MIN) then
            if (single_server_values(k) .lt. reduction_progress%values(k)) &
                 reduction_progress%values(k)=single_server_values(k)
          else if (reduction_op == MAX) then
            if (single_server_values(k) .gt. reduction_progress%values(k)) &
                 reduction_progress%values(k)=single_server_values(k)
          end if
        end do
      end if
    end if
  end subroutine integrate_io_server_collective_values

  !> Handles the case where the local MONC processes have completed their collective operation for a specific reduction
  !! and, for this IO server, it either needs to send its value to the master IO server or, if it is the master, check
  !! for completion
  !! @param io_configuration Configuration state of the IO server
  !! @param reduction_progress The specific reduction progress data item that represents this reduction
  !! @param z_size Size in Z
  !! @param reduction_progress_location Location in the reduction progresses list that this single progress item resides at
  subroutine handle_local_moncs_completed_collective(io_configuration, has_root, reduction_progress)
    type(io_configuration_type), intent(inout) :: io_configuration
    logical, intent(in) :: has_root
    type(reduction_progress_type), pointer, intent(inout) :: reduction_progress

    integer :: target_rank

    target_rank = merge(reduction_progress%root, EDAT_ALL, has_root)
    call edatFireEvent(reduction_progress%values, EDAT_DOUBLE, size(reduction_progress%values), &
          target_rank, trim(generate_reduction_key(reduction_progress%field_name, &
          reduction_progress%timestep, reduction_progress%reduction_operator)))

    if (has_root .and. io_configuration%my_io_rank .ne. reduction_progress%root) then
      ! Clean up the state here as don't need it further
      call check_thread_status(forthread_rwlock_wrlock(reduction_progress_rwlock))
      call c_remove(reduction_progresses, generate_reduction_key(reduction_progress%field_name, &
        reduction_progress%timestep, reduction_progress%reduction_operator))
      call check_thread_status(forthread_rwlock_unlock(reduction_progress_rwlock))
      call check_thread_status(forthread_mutex_destroy(reduction_progress%mutex))
      if (allocated(reduction_progress%values)) deallocate(reduction_progress%values)
      deallocate(reduction_progress)
    end if
  end subroutine handle_local_moncs_completed_collective

  subroutine reduction_procedure(events, number_events)
	  type(c_ptr), intent(in), target :: events
    integer(c_int), value, intent(in) :: number_events

    type(reduction_progress_type), pointer :: reduction_progress
    character(len=STRING_LENGTH) :: field_name
    integer :: timestep, reduction_op, i
    type(EDAT_Event) :: processed_events(number_events)

    call getEvents(events, number_events,  processed_events)

    call unpackage_reduction_key(processed_events(1)%metadata%event_id, field_name, timestep, reduction_op)
    reduction_progress=>find_reduction_progress(timestep, reduction_op, trim(field_name))

    do i=1, number_events
      if (processed_events(i)%metadata%source .ne. edatGetRank()) then
        call integrate_io_server_collective_values(reduction_op, reduction_progress, &
          processed_events(i)%double_data, processed_events(i)%metadata%number_elements, .false.)
      end if
    end do
    if (reduction_progress%reduction_operator == MEAN) then
      reduction_progress%values=reduction_progress%values/edatGetNumRanks()
    end if

    call reduction_progress%completion_procedure(this_io_configuration, reduction_progress%values, &
      field_name, timestep)

    call check_thread_status(forthread_rwlock_wrlock(reduction_progress_rwlock))
    call c_remove(reduction_progresses, generate_reduction_key(field_name, timestep, reduction_op))
    call check_thread_status(forthread_rwlock_unlock(reduction_progress_rwlock))
    call check_thread_status(forthread_mutex_destroy(reduction_progress%mutex))
    if (allocated(reduction_progress%values)) deallocate(reduction_progress%values)
    deallocate(reduction_progress)
  end subroutine reduction_procedure

  !> Finds or adds a specific reduction progress based upon the timestep and reduction operator. If none can be found
  !! then a new progress is added in. With new progresses this procedure will initialise them
  !! @param timestep The timestep to match
  !! @param reduction_operator The reduction operator to match
  !! @param field_name The name of the field that the reduction type represents
  !! @param num_vectors The number of reduction vectors (items to reduce) to be stored
  !! @returns A reduction progress data object
  function find_or_add_reduction_progress(timestep, reduction_operator, has_root, root, &
      my_rank, field_name, completion_procedure)
    integer, intent(in) :: timestep, reduction_operator, root, my_rank
    logical, intent(in) :: has_root
    type(reduction_progress_type), pointer :: find_or_add_reduction_progress
    character(len=*), intent(in) :: field_name
    procedure(handle_completion), optional :: completion_procedure

    class(*), pointer :: generic
    type(reduction_progress_type), pointer :: new_progress

    find_or_add_reduction_progress=>find_reduction_progress(timestep, reduction_operator, field_name)
    if (.not. associated(find_or_add_reduction_progress)) then
      call check_thread_status(forthread_rwlock_wrlock(reduction_progress_rwlock))
      find_or_add_reduction_progress=>find_reduction_progress(timestep, reduction_operator, field_name, issue_read_lock=.false.)
      if (.not. associated(find_or_add_reduction_progress)) then
        allocate(new_progress)
        call check_thread_status(forthread_mutex_init(new_progress%mutex, -1))
        new_progress%timestep=timestep
        new_progress%reduction_operator=reduction_operator
        new_progress%contributed_moncs=0
        new_progress%root=root
        new_progress%field_name=field_name

        if (root == my_rank .or. .not. has_root) then
          call edatSubmitTask(reduction_procedure, 1, EDAT_ALL, &
            trim(generate_reduction_key(field_name, timestep, reduction_operator)))
        end if

        if (present(completion_procedure)) then
          new_progress%completion_procedure=>completion_procedure
        else
          new_progress%completion_procedure=>null()
        end if
        generic=>new_progress
        call c_put_generic(reduction_progresses, generate_reduction_key(field_name, timestep, reduction_operator), &
             generic, .false.)
        find_or_add_reduction_progress=>new_progress
      end if
      call check_thread_status(forthread_rwlock_unlock(reduction_progress_rwlock))
    end if
    if (.not. associated(find_or_add_reduction_progress%completion_procedure) .and. present(completion_procedure)) then
      find_or_add_reduction_progress%completion_procedure=>completion_procedure
    end if
  end function find_or_add_reduction_progress

  !> Locates a specific reduction progress based upon the timestep, operator and field name
  !! @param timestep The timestep to search for
  !! @param reduction_operator The reduction operator to search for
  !! @param field_name The field name which must match
  !! @param reduction_progress_location Optional location which is set to be the index of the matching progress item
  !! @returns Pointer to the reduction progress or null if none is found
  function find_reduction_progress(timestep, reduction_operator, field_name, issue_read_lock)
    integer, intent(in) :: timestep, reduction_operator
    logical, intent(in), optional :: issue_read_lock
    type(reduction_progress_type), pointer :: find_reduction_progress
    character(len=*), intent(in) :: field_name

    class(*), pointer :: generic
    logical :: do_read_lock

    if (present(issue_read_lock)) then
      do_read_lock=issue_read_lock
    else
      do_read_lock=.true.
    end if

    if (do_read_lock) call check_thread_status(forthread_rwlock_rdlock(reduction_progress_rwlock))
    generic=>c_get_generic(reduction_progresses, generate_reduction_key(field_name, timestep, reduction_operator))
    if (do_read_lock) call check_thread_status(forthread_rwlock_unlock(reduction_progress_rwlock))
    if (associated(generic)) then
      select type(generic)
      type is (reduction_progress_type)
        find_reduction_progress=>generic
      end select
    else
      find_reduction_progress=>null()
    end if
  end function find_reduction_progress

  !> Generates the lookup key that is used for the map storage of reduction progresses
  !! @param field_name The field name
  !! @param timestep The timestep
  !! @param reduction_operator The reduction operator
  character(len=STRING_LENGTH) function generate_reduction_key(field_name, timestep, reduction_operator)
    character(len=*), intent(in) :: field_name
    integer, intent(in) :: timestep, reduction_operator

    generate_reduction_key=trim(field_name)//"#"//trim(conv_to_string(timestep))//"#"// trim(conv_to_string(reduction_operator))
  end function generate_reduction_key


  subroutine unpackage_reduction_key(key, field_name, timestep, reduction_op)
    character(len=STRING_LENGTH), intent(in) :: key
    character(len=STRING_LENGTH), intent(out) :: field_name
    integer, intent(out) :: timestep, reduction_op

    integer :: first_hash, second_hash

    first_hash=scan(key, "#")
    second_hash=scan(key(first_hash+1:), "#") + first_hash

    field_name=key(1:first_hash-1)
    timestep=conv_to_integer(key(first_hash+1:second_hash-1))
    reduction_op=conv_to_integer(key(second_hash+1:))
  end subroutine unpackage_reduction_key

  !> Given the map of action attributes this procedure will identify the reduction operator that has been
  !! selected by the configuration
  !! @param action_attributes Action attributes from the IO server configuration
  !! @returns The reduction operator
  integer function get_reduction_operator(op_string)
    character(len=*), intent(in) :: op_string

    if (op_string .eq. "mean") then
      get_reduction_operator=MEAN
    else if (op_string .eq. "min") then
      get_reduction_operator=MIN
    else if (op_string .eq. "max") then
      get_reduction_operator=MAX
    else if (op_string .eq. "sum") then
      get_reduction_operator=SUM
    else
      call log_log(LOG_ERROR, "Reduction operator '"//trim(op_string)//"' not recognised")
    end if
  end function get_reduction_operator
end module reduction_inter_io_mod
