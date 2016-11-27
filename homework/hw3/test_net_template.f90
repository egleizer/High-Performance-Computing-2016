program test_network
	    use network
	    implicit none
	integer :: N01,L1,Nt1,qmax
	N01=5.d0
	Nt1=8.d0
	L1=3.d0
	call generate(N01,L1,Nt1,qmax)
	call adjacency_matrix()
!subroutine calls below should be modified as required and variables should be declared above as needed.
!call generate()
!call adjacency_matrix()
!call output_anet()
!call output_enet()
!call output_q()
!call connectivity()
!call vary_connectivity()
end program test_network




