!> Abstraction layer around MPI, this issues and marshals the lower level communication details
module mpi_communication_mod
  use datadefn_mod, only : STRING_LENGTH
  use collections_mod, only : map_type, c_put_integer
  use configuration_parser_mod, only : io_configuration_data_definition_type, io_configuration_inter_communication_description
  use logging_mod, only : LOG_ERROR, log_log, log_master_log
  use io_server_client_mod, only : ARRAY_FIELD_TYPE, MAP_FIELD_TYPE, INTEGER_DATA_TYPE, BOOLEAN_DATA_TYPE, STRING_DATA_TYPE, &
       FLOAT_DATA_TYPE, DOUBLE_DATA_TYPE, COMMAND_TAG, DATA_TAG, INTER_IO_COMMUNICATION, get_data_description_from_name, &
       data_sizing_description_type, populate_type_extents
  use forthread_mod, only : forthread_mutex_lock, forthread_mutex_unlock, forthread_mutex_init, forthread_mutex_destroy
  use data_utils_mod, only : check_thread_status
  use mpi, only : MPI_COMM_WORLD, MPI_SOURCE, MPI_INT, MPI_BYTE, MPI_STATUS_SIZE, MPI_REQUEST_NULL, &
       MPI_STATUS_IGNORE, MPI_STATUSES_IGNORE, MPI_ANY_SOURCE, MPI_THREAD_MULTIPLE, MPI_THREAD_SERIALIZED
  use iso_c_binding
  use edat, only : edatLockComms, edatUnlockComms
  implicit none

#ifndef TEST_MODE
  private
#endif

  integer, parameter :: MS_WAIT_BETWEEN_TESTS=100

  !< Interface to the C usleep Linux call which allows us to sleep for a specific number of MS
  interface
     subroutine usleep(useconds) bind(C)
       use iso_c_binding
       implicit none
       integer(c_int32_t), value :: useconds
     end subroutine usleep
  end interface

  integer :: mpi_threading_mode
  integer, volatile :: mpi_mutex
  logical :: manage_mpi_thread_safety

  public free_mpi_type, get_number_io_servers, get_my_io_rank, initialise_mpi_communication, &
       lock_mpi, unlock_mpi, wait_for_mpi_request, waitall_for_mpi_requests
contains

  !> Initialises MPI communication
  !! @param provided_threading The provided threading mode
  subroutine initialise_mpi_communication(provided_threading)
    integer, intent(in) :: provided_threading

    mpi_threading_mode=provided_threading
    if (mpi_threading_mode .ne. MPI_THREAD_MULTIPLE .and. mpi_threading_mode .ne. MPI_THREAD_SERIALIZED) then
      call log_master_log(LOG_ERROR, "You must run MONC in MPI thread serialized or thread multiple mode for the IO server")
    end if
    manage_mpi_thread_safety=provided_threading == MPI_THREAD_SERIALIZED
    call check_thread_status(forthread_mutex_init(mpi_mutex, -1))
  end subroutine initialise_mpi_communication

  !> If we are explicitly managing MPI thread safety (SERIALIZED mode) then locks MPI
  subroutine lock_mpi()
    if (manage_mpi_thread_safety) call check_thread_status(forthread_mutex_lock(mpi_mutex))
    call edatLockComms()
  end subroutine lock_mpi

  !> If we are explicitly managing MPI thread safety (SERIALIZED mode) then unlocks MPI
  subroutine unlock_mpi()
    if (manage_mpi_thread_safety) call check_thread_status(forthread_mutex_unlock(mpi_mutex))
    call edatUnlockComms()
  end subroutine unlock_mpi

  !> Pauses for a specific number of ms to allow for MPI interleaving, this is to avoid starvation
  subroutine pause_for_mpi_interleaving()
    call usleep(int(MS_WAIT_BETWEEN_TESTS, c_int32_t))
  end subroutine pause_for_mpi_interleaving

  !> Waits for a specific MPI request to complete, either by managing thread safety and interleaving or just a call to MPI
  !! if we are in multiple mode
  !! @param request The MPI request handle
  subroutine wait_for_mpi_request(request, status)
    integer, intent(inout) :: request
    integer, intent(inout), optional :: status(MPI_STATUS_SIZE)

    integer :: ierr, flag

    if (manage_mpi_thread_safety) then
      flag=0
      do while (flag .ne. 1)
        call lock_mpi()
        if (present(status)) then
          call mpi_test(request, flag, status, ierr)
        else
          call mpi_test(request, flag, MPI_STATUS_IGNORE, ierr)
        end if
        call unlock_mpi()
        if (flag .ne. 1) call pause_for_mpi_interleaving()
      end do
    else
      if (present(status)) then
        call mpi_wait(request, status, ierr)
      else
        call mpi_wait(request, MPI_STATUS_IGNORE, ierr)
      end if
    end if
  end subroutine wait_for_mpi_request

  !> Waits for all MPI requests to complete, either by managing thread safety and interleaving or just a call to MPI
  !! if we are in multiple mode
  !! @param requests The MPI request handles to wait for
  !! @param count The number of request handles to wait for
  subroutine waitall_for_mpi_requests(requests, count)
    integer, dimension(:), intent(inout) :: requests
    integer, intent(in) :: count

    integer :: ierr, flag

    if (manage_mpi_thread_safety) then
      flag=0
      do while (flag .ne. 1)
        call lock_mpi()
        call mpi_testall(count, requests, flag, MPI_STATUSES_IGNORE, ierr)
        call unlock_mpi()
        if (flag .ne. 1) call pause_for_mpi_interleaving()
      end do
    else
      call mpi_waitall(count, requests, MPI_STATUSES_IGNORE, ierr)
    end if
  end subroutine waitall_for_mpi_requests

  !> Retrieves the number of IO servers that are running in total
  !! @param io_comm The IO server communicator
  !! @returns The number of running IO servers
  integer function get_number_io_servers(io_comm)
    integer, intent(in) :: io_comm

    integer :: number, ierr

    call mpi_comm_size(io_comm, number, ierr)
    get_number_io_servers=number
  end function get_number_io_servers

  !> Retrieves my IO server rank out of the number of IO servers that are running
  !! @param io_comm The IO server communicator
  !! @returns My IO server rank
  integer function get_my_io_rank(io_comm)
    integer, intent(in) :: io_comm

    integer :: number, ierr

    call mpi_comm_rank(io_comm, number, ierr)
    get_my_io_rank=number
  end function get_my_io_rank

  !> Frees an MPI type, used in clean up
  !! @param the_type The MPI type to free up
  subroutine free_mpi_type(the_type)
    integer, intent(in) :: the_type

    integer :: ierr

    call mpi_type_free(the_type, ierr)
  end subroutine free_mpi_type
end module mpi_communication_mod
