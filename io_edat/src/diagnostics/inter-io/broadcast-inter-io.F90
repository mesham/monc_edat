!> Broadcast inter IO communication which sends a value from one IO server to all others. This tracks field name and timestep
!! and only issues one call (and one results call to completion) for that combination
module broadcast_inter_io_mod
  use datadefn_mod, only : DEFAULT_PRECISION, DOUBLE_PRECISION, STRING_LENGTH
  use configuration_parser_mod, only : io_configuration_type, io_configuration_inter_communication_description
  use collections_mod, only : hashmap_type, list_type, iterator_type, mapentry_type, c_add_string, c_remove, c_free, &
       c_get_generic, c_get_string, c_put_generic, c_generic_at, c_get_iterator, c_has_next, c_next_mapentry, c_next_string, &
       c_is_empty
  use conversions_mod, only : conv_to_string, conv_to_integer
  use forthread_mod, only : forthread_mutex_init, forthread_mutex_lock, forthread_mutex_trylock, forthread_mutex_unlock, &
       forthread_mutex_destroy, forthread_rwlock_rdlock, forthread_rwlock_wrlock, forthread_rwlock_unlock, &
       forthread_rwlock_init, forthread_rwlock_destroy
  use data_utils_mod, only : check_thread_status

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

  !< Type keeping track of broadcast statuses
  type inter_io_broadcast
     logical :: handled
     integer :: mutex
     procedure(handle_completion), pointer, nopass :: completion_procedure
  end type inter_io_broadcast

  type(io_configuration_type), pointer :: this_io_configuration
  type(hashmap_type), volatile :: broadcast_statuses
  integer, volatile :: broadcast_statuses_rwlock
  logical, volatile :: initialised=.false.

  public init_broadcast_inter_io, perform_inter_io_broadcast, finalise_broadcast_inter_io
contains

  !> Initialises the broadcast inter IO functionality
  !! @param io_configuration The IO server configuration
  subroutine init_broadcast_inter_io(io_configuration)
    type(io_configuration_type), intent(inout), target :: io_configuration

    if (.not. initialised) then
      this_io_configuration=>io_configuration
      initialised=.true.
      call check_thread_status(forthread_rwlock_init(broadcast_statuses_rwlock, -1))
    end if
  end subroutine init_broadcast_inter_io

  !> Finalises the broadcast inter IO functionality
  subroutine finalise_broadcast_inter_io()

    if (initialised) then
      call check_thread_status(forthread_rwlock_destroy(broadcast_statuses_rwlock))
      initialised=.false.
    end if
  end subroutine finalise_broadcast_inter_io

  !> Performs an inter IO broadcast of data from the root to all other IO servers. Note that this is on the IO server (and not
  !! MONC level) so might require some translation between the user's logical view and this view. Broadcasts are only issued once
  !! for a specific field_name and timestep pair.
  !! @param io_configuration Configuration of the IO server
  !! @param field_values The values to communicate
  !! @param field_size Number of elements to communicate
  !! @param field_name Field name that the reduction will be performed over
  !! @param root The root IO server process
  !! @param timestep The timestep this is issued at
  !! @param completion_procedure Callback completion procedure
  subroutine perform_inter_io_broadcast(io_configuration, field_values, field_size, field_name, root, &
       timestep, completion_procedure)
    type(io_configuration_type), intent(inout) :: io_configuration
    real(kind=DOUBLE_PRECISION), dimension(:) :: field_values
    integer, intent(in) :: field_size, root, timestep
    character(len=*), intent(in) :: field_name
    procedure(handle_completion) :: completion_procedure

    type(inter_io_broadcast), pointer :: broadcast_item

    broadcast_item=>find_or_add_broadcast_item(field_name, timestep, completion_procedure)

    call check_thread_status(forthread_mutex_lock(broadcast_item%mutex))
    if (.not. broadcast_item%handled) then
      broadcast_item%handled=.true.
      call check_thread_status(forthread_mutex_unlock(broadcast_item%mutex))
      call edatSubmitTask(broadcast_procedure, 1, root, trim(field_name)//"#"//trim(conv_to_string(timestep)))

      if (io_configuration%my_io_rank == root) then
        call edatFireEvent(field_values, EDAT_DOUBLE, size(field_values), &
          EDAT_ALL, trim(field_name)//"#"//trim(conv_to_string(timestep)))
      end if
    else
      call check_thread_status(forthread_mutex_unlock(broadcast_item%mutex))
    end if
  end subroutine perform_inter_io_broadcast

  subroutine broadcast_procedure(events, number_events)
	  type(c_ptr), intent(in), target :: events
    integer(c_int), value, intent(in) :: number_events

    character(len=STRING_LENGTH) :: field_name
    integer :: timestep, hash_loc
    type(EDAT_Event) :: processed_events(number_events)
    type(inter_io_broadcast), pointer :: broadcast_item

    call getEvents(events, number_events,  processed_events)

    hash_loc=scan(processed_events(1)%metadata%event_id, "#")
    field_name=processed_events(1)%metadata%event_id(:hash_loc-1)
    timestep=conv_to_integer(trim(processed_events(1)%metadata%event_id(hash_loc+1:)))

    broadcast_item=>find_broadcast_item(field_name, timestep, .true.)

    if (associated(broadcast_item%completion_procedure)) then
      call broadcast_item%completion_procedure(this_io_configuration, processed_events(1)%double_data, trim(field_name), timestep)
    end if
    call check_thread_status(forthread_rwlock_wrlock(broadcast_statuses_rwlock))
    call c_remove(broadcast_statuses, trim(processed_events(1)%metadata%event_id))
    call check_thread_status(forthread_rwlock_unlock(broadcast_statuses_rwlock))
    call check_thread_status(forthread_mutex_destroy(broadcast_item%mutex))
    deallocate(broadcast_item)
  end subroutine broadcast_procedure

  !> Locates and returns or adds and returns a specific broadcast item representing a timestep and field
  !! @param field_name The field name this represents
  !! @param timestep The timestep this represents
  !! @param completion_procedure The (optional) completion procedure which is called once values are received
  !! @returns The existing or new broadcast item
  function find_or_add_broadcast_item(field_name, timestep, completion_procedure)
    character(len=*), intent(in) :: field_name
    integer, intent(in) :: timestep
    procedure(handle_completion), optional :: completion_procedure
    type(inter_io_broadcast), pointer :: find_or_add_broadcast_item

    class(*), pointer :: generic

    find_or_add_broadcast_item=>find_broadcast_item(field_name, timestep, .true.)
    if (.not. associated(find_or_add_broadcast_item)) then
      call check_thread_status(forthread_rwlock_wrlock(broadcast_statuses_rwlock))
      find_or_add_broadcast_item=>find_broadcast_item(field_name, timestep, .false.)
      if (.not. associated(find_or_add_broadcast_item)) then
        allocate(find_or_add_broadcast_item)
        if (present(completion_procedure)) then
          find_or_add_broadcast_item%completion_procedure=>completion_procedure
        else
          find_or_add_broadcast_item%completion_procedure=>null()
        end if
        find_or_add_broadcast_item%handled=.false.
        call check_thread_status(forthread_mutex_init(find_or_add_broadcast_item%mutex, -1))
        generic=>find_or_add_broadcast_item
        call c_put_generic(broadcast_statuses, trim(field_name)//"#"//conv_to_string(timestep), generic, .false.)
      end if
      call check_thread_status(forthread_rwlock_unlock(broadcast_statuses_rwlock))
    end if
  end function find_or_add_broadcast_item

  !> Finds a specific broadcast item or null if none is found
  !! @param field_name Corresponding field name to find
  !! @param timestep Corresponding timestep to find
  !! @param do_read_lock Whether to issue a read lock or not
  !! @returns The corresponding broadcast status item or null if none is found
  function find_broadcast_item(field_name, timestep, do_read_lock)
    character(len=*), intent(in) :: field_name
    integer, intent(in) :: timestep
    logical, intent(in) :: do_read_lock
    type(inter_io_broadcast), pointer :: find_broadcast_item

    class(*), pointer :: generic

    if (do_read_lock) call check_thread_status(forthread_rwlock_rdlock(broadcast_statuses_rwlock))
    generic=>c_get_generic(broadcast_statuses, trim(field_name)//"#"//conv_to_string(timestep))
    if (do_read_lock) call check_thread_status(forthread_rwlock_unlock(broadcast_statuses_rwlock))

    if (associated(generic)) then
      select type(generic)
        type is (inter_io_broadcast)
          find_broadcast_item=>generic
      end select
    else
      find_broadcast_item=>null()
    end if
  end function find_broadcast_item
end module broadcast_inter_io_mod
