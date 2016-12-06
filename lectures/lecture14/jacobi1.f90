!Solve 1-d heat equation with Jacobi iteration
program jacobi1
	implicit none
	integer :: i1,k1,n,kmax
	real(kind=8) :: a,b !b.c.'s
    real(kind=8) :: dx,dx2f,error !grid spacing
	real(kind=8), allocatable, dimension(:) :: x,S,T,Tnew !grid, source, Temperatures


    open(unit=12,file='data.in')
    read(12,*) n
    read(12,*) a
    read(12,*) b


    allocate(x(0:n),S(0:n),T(0:n),Tnew(0:n))

    kmax = 1000

    !grid
    dx = 1.d0/dble(n+1)
    dx2f = 0.5*dx**2

    do i1=0,n+1
        x(i1) = i1*dx
    end do

    !set initial condition
    T = (b-a)*x + a

    !set source function
    S = exp(-100.d0*(x-0.5)**2)

    do k1=1,kmax
        Tnew(1:n) = (-S(1:n) - (T(0:n-1) + T(2:n+1)))*dx2f
        error = maxval(abs(Tnew(1:n)-T(1:n)))
        print *, 'k,error=',k1,error
        T(1:n)=Tnew(1:n)
    end do


end program jacobi1

