!module for computing circumference and area of circle
module circle
    implicit none
    real(kind=8) :: pi
    save
contains

subroutine initialize_pi()
    implicit none

    pi = acos(-1.d0)

end subroutine initialize_pi

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

end module circle




