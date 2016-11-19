!Illustrates use of parallel region

program firstomp
	use omp_lib !makes OpenMP routines, variables available
	implicit none
	integer :: NumThreads,threadID


!$OMP PARALLEL PRIVATE(threadID)
	NumThreads = omp_get_num_threads()
	threadID = omp_get_thread_num()

    if (threadID==0) then
        call subroutine1(in1,out1)
    elseif (threadID==1) then
        call subroutine1(in2,out2)
    end if

!$OMP END PARALLEL

    print *, 'exited parallel region, this is thread', threadID

end program firstomp




