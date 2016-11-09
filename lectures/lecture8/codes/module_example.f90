!illustration of module structure
module module_example2

    !1. Variable declarations, these variables are available 
    !to all function and subroutines contained in module
    implicit none
    integer :: i1,j1
    save

contains
!2. Functions and subroutines

real(kind=8) function f1(x)
    !compute product of x and i1
    !assumes i1 has been initialized elsewhere
    real(kind=8) x

    f1 = x*i1
end function f1

subroutine f2(x,y)
    !scales input variables with j1
    !assumes j1 has been initialized elsewhere
    real(kind=8) :: x,y

    x = x*j1
    y = y*j1

end subroutine f2

end module module_example2




module module_name

    !1. variable declarations

contains

    !2. subroutines and functions

end module module_name


program module_example
    use module_name
    implicit none
    !variable decilarations

    !code

end program module_example





