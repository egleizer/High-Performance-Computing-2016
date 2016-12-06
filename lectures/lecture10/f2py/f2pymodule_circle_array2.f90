!module for computing circumference, area, and "mass" of circle
module circle
    implicit none
    real(kind=8) :: pi
    real(kind=8), allocatable, dimension(:) :: weights
    save
contains

subroutine initialize_pi()
    implicit none

    pi = acos(-1.d0)

end subroutine initialize_pi

subroutine initialize_weights(N)
    implicit none
    integer, intent(in) :: N
    allocate(weights(N))
    weights = 1.d0
end subroutine initialize_weights


real(kind=8) function circumference(radius)
    !compute circumference of circle given the radius
    implicit none
    real(kind=8), intent(in) :: radius

    circumference = 2.d0*pi*radius

end function circumference


real(kind=8) function area(radius)
    !compute area of circle given the radius
    implicit none
    real(kind=8), intent(in) :: radius

    area = pi*(radius**2)

end function area

subroutine compute_mass(N,radius,mass)
    !compute mass = weights*area
    implicit none
    integer, intent(in) :: N
    real(kind=8), intent(in) :: radius
    real(kind=8), intent(out),dimension(N) :: mass

    if (allocated(weights)) then
        mass = weights*area(radius)
    else
        allocate(weights(1))
        weights = 1.d0
        print *, 'error, weights not allocated'
    end if
end subroutine compute_mass

end module circle




