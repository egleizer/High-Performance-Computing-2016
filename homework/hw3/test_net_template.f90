program test_network
<<<<<<< HEAD
	    use network
	    implicit none
	integer :: N01,L1,Nt1,qmax
	N01=5.d0
	Nt1=8.d0
	L1=3.d0
	call generate(N01,L1,Nt1,qmax)
	call adjacency_matrix()
=======
	use network
	implicit none


>>>>>>> 12183415aaa9338c10fc8c6a264bc3c5ba57c88a
!subroutine calls below should be modified as required and variables should be declared above as needed.
!call generate()
!call adjacency_matrix()
!call output_anet()
!call output_enet()
!call output_q()
!call connectivity()
!call vary_connectivity()
<<<<<<< HEAD
end program test_network


=======


end program test_network


subroutine output_anet()
	!creates data file anet.dat which can be read in python 
	!with f = np.loadtxt("anet.dat")
	!Works for 10000 x 10000 matrices and smaller
	use network
	implicit none
	integer :: i1,n
	n = size(anet,1)
	open(unit=11,file='anet.dat')
	do i1 = 1,n 
		write(11,'(10000I4)') anet(i1,:)
	end do
	close(11)
end subroutine output_anet()


subroutine output_enet()
	!creates data file enet.dat which can be read in python 
	!with f = np.loadtxt("enet.dat")
	use network
	implicit none
	integer :: i1,n
	n = size(enet,1)
	open(unit=11,file='enet.dat')
	do i1 = 1,n 
		write(11,*) enet(i1,:)
	end do
	close(11)
end subroutine output_enet()

subroutine output_qnet()
	!creates data file qnet.dat which can be read in python 
	!with f = np.loadtxt("qnet.dat")
	use network
	implicit none
	integer :: i1,n
	n = size(qnet,1)
	open(unit=11,file='qnet.dat')
	do i1 = 1,n 
		write(11,*) qnet(i1,:)
	end do
	close(11)
end subroutine output_qnet()
>>>>>>> 12183415aaa9338c10fc8c6a264bc3c5ba57c88a


