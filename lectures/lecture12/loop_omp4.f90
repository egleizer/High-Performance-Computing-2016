!Simple example of using OpenMP with nested loops
program loop_omp3
	use omp_lib !makes OpenMP routines, variables available
	implicit none
	integer :: NumThreads,threadID
    integer :: i1,j1,N,M
    real(kind=8), dimension(4,5) :: x,y,z


    M=size(x,1)
    N=size(x,2)

!initialize y,z
    do j1=1,N
        do i1=1,M
            y(i1,j1)=i1*j1
        end do
    end do
    z = sin(y)
    x(:,1) = z(:,1)

!$OMP parallel do private(j1)
do i1 = 1,M
    do j1 = 2,N
        x(i1,j1) = x(i1,j1-1)
    end do
end do
!$OMP end parallel do

!print *, 'test:', maxval(abs(x-y-z))
print *, 'x=',x
end program loop_omp3




