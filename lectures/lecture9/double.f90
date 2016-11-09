!program and subroutine illustrating fortran
!treatment of subroutine input/output.
!Two numbers are sent into a subroutine and
!their values are doubled.

program main
	implicit none
	real(kind=8) :: x,y


	x=2.d0
	y=3.d0

	print *, 'before subroutine, x,y=',x,y
	
	call double_sub(x,y)

	print *, 'after subroutine, x,y=',x,y


end program main

subroutine double_sub(a,b)
	implicit none
	real(kind=8), intent(inout) :: a,b

	a = 2.d0*a
	b = 2.d0*b

end subroutine double_sub







