!Template code for Homework 3
module network




subroutine generate(N0,L,Nt,qmax)
	!Generate recursive matrix corresponding
	!to model parameters provided as input
	implicit none
	integer, intent(in) :: N0,L,Nt
	integer, intent(out) :: qmax

end subroutine generate


subroutine adjacency_matrix()
	!Create adjacency matrix corresponding 
	!to pre-existing edge list 


end subroutine adjacency_matrix

subroutine connectivity(c)
	!Compute connectivity of pre-computed network
	implicit none
	real(kind=8), intent(out) :: c

end subroutine connectivity


subroutine vary_connectivity(N0,Lmax,Nt,carray)
	!Compute connectivity for networks with L varying
	!between 1 and Lmax
	implicit none
	integer, intent(in) :: N0,Lmax,Nt
	real(kind=8), dimension(Lmax),intent(out) :: carray


end subroutine vary_connectivity()


end module network
