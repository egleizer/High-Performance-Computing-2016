!Getting started with OpenMP

program firstomp
	use omp_lib !makes OpenMP routines, variables available
	implicit none
	integer :: NumThreads,threadID


    NumThreads = omp_get_num_threads()
    print *, 'before parallel region, NumThreads=', NumThreads


!$OMP PARALLEL PRIVATE(threadID)
	NumThreads = omp_get_num_threads()
	threadID = omp_get_thread_num()
	print *, 'this is thread',threadID, ' of ', NumThreads
!$OMP END PARALLEL


	print *, 'exited parallel region, this is thread', threadID

end program firstomp




