!Lab 5 solution code
!Compute the integral of 4/(1+x^2) from 0 to 1 using the midpoint or trapezoid rule
!Input: number of intervals, N
!Output (to screen): number of intervals, estimated value of integral, and error
!To compile this code: gfortran -o quad.exe quadsoln.f90
!To run: ./quad.exe
!For f2py: f2py -c quadsoln.f90 -m q
!which generates q.so

program quad
	implicit none
    integer :: i1,N,quad_type
    real(kind=8) :: pi, I, error
    real(kind=8), allocatable, dimension(:) :: x,f !Task 4


    !read data from data.in
	open(unit=10, file='data.in')
        read(10,*) N
        read(10,*) quad_type
	close(10)

    !compute integral
    if (quad_type==1) then
        call midpoint(N,I)
    else
        call trapezoid(N,I)
    end if

    pi = asin(1.d0)*2.d0
    error = abs(I - pi)
    print *, 'N=', N
    print *, 'integral=', I
    print *, 'error=', error

    !Task 4: Complete code below (and add appropriate variable declarations above)
    !Write x,f(x) to file 
    allocate(x(N+1),f(N+1))
    do i1=1,N+1
        x(i1) = dble(i1-1)/dble(N)
    end do
    f = 4.d0/(1.d0 + x*x)
    open(unit=11,file='integrand.dat')
    do i1 = 1,N+1
        write(11,*) x(i1),f(i1)
    end do
    close(11)
    deallocate(x,f)
    ! To load in python:
    ! F = np.loadtxt("integrand.dat")


end program quad

!--------------------------------------------------------------------

!----------------------------------
!subroutine midpoint
! Use midpoint rule to compute 
! integral of, 4.0/(1+x^2)
! from x=0 to x=1
!----------------------------------
subroutine midpoint(N,I)
    implicit none
    integer, intent(in) :: N
    real(kind=8), intent(out) :: I
    integer :: i1
    real(kind=8) :: dx, sum_i, xm, f, a, error

     I = 0.d0 !initialize integral
     dx = 1.d0/dble(N) !interval size

    !loop over intervals computing each interval's contribution to integral
    do i1 = 1,N
        xm = dx*(dble(i1)-0.5d0) !midpoint of interval i1
        call integrand(xm,f)
        sum_i = dx*f
        I = I + sum_i !add contribution from interval to total integral
    end do

end subroutine midpoint
!--------------------------------------------------------------------

!----------------------------------
!subroutine trapezoid
! Use trapezoid rule to compute 
! integral of, 4.0/(1+x^2)
! from x=0 to x=1
!----------------------------------
subroutine trapezoid(N,I)
    implicit none
    integer, intent(in) :: N
    real(kind=8), intent(out) :: I
    integer :: i1
    real(kind=8) :: dx, sum_i, xm, f, a, error

     I = 0.d0 !initialize integral
     dx = 1.d0/dble(N) !interval size

    !loop over intervals computing each interval's contribution to integral
    do i1 = 1,N+1 
        xm = dx*(dble(i1)-1.d0) 
        call integrand(xm,f)
        if ((i1==1) .or. (i1==N)) f=f/2.d0 !adjust endpoints
        sum_i = dx*f
        I = I + sum_i !add contribution from interval to total integral
    end do

end subroutine trapezoid


!--------------------------------------------------------------------

!----------------------------------
!subroutine integrand
!   compute integrand, 4.0/(1+a^2)
!----------------------------------

subroutine integrand(a,f)
    implicit none
    real(kind=8), intent(in) :: a
    real(kind=8), intent(out) :: f
    f = 4.d0/(1.d0 + a*a)
end subroutine integrand

!--------------------------------------------------------------------



