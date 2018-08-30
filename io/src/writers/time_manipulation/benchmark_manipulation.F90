module benchmark_time_manipulation_mod
  use datadefn_mod, only : DEFAULT_PRECISION, STRING_LENGTH
  use collections_mod, only : hashmap_type, c_put_real, c_get_real, c_contains
  use conversions_mod, only : conv_single_real_to_double, conv_to_integer
  use forthread_mod, only : forthread_mutex_init, forthread_mutex_lock, forthread_mutex_unlock, forthread_mutex_destroy
  use configuration_parser_mod, only : data_values_type
  use forthread_mod, only : forthread_mutex_init, forthread_mutex_lock, forthread_mutex_unlock
  use threadpool_mod, only : check_thread_status
	use mpi_communication_mod, only : lock_mpi, unlock_mpi
  use mpi
  implicit none

#ifndef TEST_MODE
  private
#endif

	integer, volatile :: tot_vals=0, mtx, latency_count=0, io_rank, communicator
  double precision, volatile :: start_time, end_time, latency_avg=0.0, latency_min=999.9, latency_max=0.0
  logical, volatile :: started = .false.

  public is_benchmark_time_manipulation_ready_to_write, perform_benchmark_time_manipulation, &
    finalise_benchmark_time_manipulation, initialise_benchmark_time_manipulation
contains

   subroutine initialise_benchmark_time_manipulation(myio_rank, comm)
		integer, intent(in) :: myio_rank, comm

		io_rank=myio_rank
		communicator=comm

		call check_thread_status(forthread_mutex_init(mtx, -1))
   end subroutine initialise_benchmark_time_manipulation

  subroutine finalise_benchmark_time_manipulation()
	
	integer :: global_vals, global_latency_count, ierr
    double precision :: timediff, g_timediff, global_latency_avg, global_latency_min, global_latency_max

    if (.not. started) then
      timediff=0.0
    else
      timediff=end_time-start_time
    end if

    call lock_mpi()

    call MPI_Reduce(tot_vals, global_vals, 1, MPI_INTEGER, MPI_SUM, 0, communicator, ierr)
    call MPI_Reduce(timediff, g_timediff, 1, MPI_DOUBLE, MPI_MAX, 0, communicator, ierr)
    
    call MPI_Reduce(latency_avg, global_latency_avg, 1, MPI_DOUBLE, MPI_SUM, 0, communicator, ierr)
    call MPI_Reduce(latency_min, global_latency_min, 1, MPI_DOUBLE, MPI_MIN, 0, communicator, ierr)
    call MPI_Reduce(latency_max, global_latency_max, 1, MPI_DOUBLE, MPI_MAX, 0, communicator, ierr)
    call MPI_Reduce(latency_count, global_latency_count, 1, MPI_INTEGER, MPI_SUM, 0, communicator, ierr)

    call unlock_mpi()

		if (io_rank == 0) then
    	print *, "Bandwidth: ", MPI_Wtime() - start_time, tot_vals
    	print *, "Average Latency: ", latency_avg/latency_count, " based on ", latency_count, " items"
			print *, "Minimum Latency: ", latency_min
			print *, "Maximum Latency: ", latency_max
		end if
  end subroutine finalise_benchmark_time_manipulation

  logical function is_benchmark_time_manipulation_ready_to_write(latest_time, output_frequency, write_time, &
       latest_timestep, write_timestep)
    real, intent(in) :: latest_time, output_frequency, write_time
    integer, intent(in) :: latest_timestep, write_timestep

    is_benchmark_time_manipulation_ready_to_write=latest_timestep .ge. write_timestep
  end function is_benchmark_time_manipulation_ready_to_write

  type(data_values_type) function perform_benchmark_time_manipulation(instant_values, output_frequency, &
       field_name, timestep, time)
    real(kind=DEFAULT_PRECISION), dimension(:), intent(in) :: instant_values
    real, intent(in) :: output_frequency
    real(kind=DEFAULT_PRECISION), intent(in) :: time
    character(len=*), intent(in) :: field_name
    integer, intent(in) :: timestep

		double precision :: time_diff

    call check_thread_status(forthread_mutex_lock(mtx))

    if (index(field_name, "timestamp") .gt. 0) then
			time_diff=MPI_Wtime()-instant_values(1)
      latency_avg=latency_avg+time_diff
      latency_count=latency_count+1
      if (latency_min .gt. time_diff) latency_min=time_diff
			if (latency_max .lt. time_diff) latency_max=time_diff
    end if

    if (.not. started) then
      started=.true.
      start_time=MPI_Wtime()
    end if
    tot_vals=tot_vals+1
		end_time=MPI_Wtime()
    call check_thread_status(forthread_mutex_unlock(mtx))
  end function perform_benchmark_time_manipulation
end module benchmark_time_manipulation_mod
