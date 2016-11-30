!module to use RK4 time marching to solve an initial value problem
!Solves: dy/dt = RHS(y,t)
!To compile: gfortran -o ode_test.exe ode.f90
!To run: ./ode_test.exe
module ode

contains
!---------------------------
subroutine rk4(t0,y0,dt,nt,y)
    !4th order RK method
    !input:
    !t0: initial time
    !y0: initial condition (array)
    !dt: time step
    !nt: number of time steps
    !output:
    !y: solution at time t+dt*nt (array)
    implicit none
    real(kind=8), dimension(:), intent(in) :: y0
    real(kind=8), intent(in) :: t0,dt
	integer, intent (in) :: nt
    real(kind=8), dimension(size(y0)), intent(out) :: y
    real(kind=8), dimension(size(y0)) :: f1, f2, f3, f4
    real(kind=8) :: t,halfdt,fac
	integer:: k

        halfdt = 0.5d0*dt
        fac = 1.d0/6.d0

        y = y0 !initial condition
        t = t0 !initial time

        do k = 1, nt !advance nt time steps

           f1 = dt*RHS(t, y)

           f2 = dt*RHS(t + halfdt, y + 0.5d0*f1)

           f3 = dt*RHS(t + halfdt, y + 0.5d0*f2)

           f4 = dt*RHS(t + dt, y + f3)

           y = y + (f1 + 2*f2  + 2*f3 + f4)*fac

           t = t + dt*dble(k)

        end do
end subroutine rk4
!---------------------------
subroutine euler(t0,y0,dt,nt,y)
    !explicit Euler method
    !input:
    !t0: initial time
    !y0: initial condition (array)
    !dt: time step
    !nt: number of time steps
    !output:
    !y: solution at time t+dt*nt (array)
    implicit none
    real(kind=8), dimension(:), intent(in) :: y0
    real(kind=8), intent(in) :: t0,dt
	integer, intent (in) :: nt
    real(kind=8), dimension(size(y0)), intent(out) :: y
    real(kind=8) :: t,halfdt,fac
	integer:: k

    y = y0 !initial condition
    t = t0 !initial time
    do k = 1,nt !advance nt time steps

        y = y + dt*RHS(t,y)
        t = t + dt

    end do


end subroutine euler
!---------------------------
function RHS(t,f)
    !RHS for dy_i/dt = -y_i
    implicit none
    real(kind=8), intent(in) :: t
    real(kind=8), dimension(:), intent(in) :: f
    real(kind=8), dimension(size(f)) :: RHS

    RHS = -f
end function RHS
!---------------------------

end module ode


program test_rk4
    !test ode solvers with dy/dt = -y, y(t=0) = sin(x)
    use ode
    implicit none
    integer :: n,i1
    real(kind=8), allocatable, dimension(:) :: x,f0,feuler,frk4



    open(unit=10,file='data.in')
    read(10,*) n
    close(10)

    allocate(x(n),f0(n),feuler(n),frk4(n))

    do i1=1,n
        x(i1) = dble(i1-1)/n
    end do

    f0 = sin(x)

    call euler(0.d0,f0,0.01d0,100,feuler)
    call rk4(0.d0,f0,0.01d0,100,frk4)

!Lab 8, part 2:  Add code to compute error:

    print *, 'euler error = ',
    print *, 'rk4 error = ',


end program test_rk4





