program test
use circle
implicit none
real(kind=8), allocatable, dimension(:) :: mass

call initialize_pi()
call initialize_weights(10)
allocate(mass(10))

print *, circumference(2.d0)
call compute_mass(10,2.d0,mass)
print *, 'mass=',mass

end program test
