!Simple example of OpenMP parallel loop
program firstomp
	use omp_lib !makes OpenMP routines, variables available
	implicit none
	integer :: NumThreads,threadID
    integer :: i1
    real(kind=8), dimension(10) :: x,y,z

    do i1=1,size(y)
        y(i1)=i1
    end do
    z = sin(y)


!$OMP PARALLEL PRIVATE(threadID)
	NumThreads = omp_get_num_threads()
	threadID = omp_get_thread_num()
	print *, 'this is thread',threadID, ' of ', NumThreads
!$OMP END PARALLEL

!$OMP parallel do
do i1 = 1,n
    x(i1) = y(i1) + z(i1)
end do
!$OMP end parallel do


    print *, 'exited parallel region, this is thread', threadID

end program firstomp




