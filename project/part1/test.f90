program test
	use network
	use rwmodule
	implicit none
	integer, dimension(405,100) :: qnetm
	integer :: qmaxm, qmax, numThreads
	real(kind=8) :: qvarm, walltime
	integer, dimension(6) :: alist1
	integer, dimension(2) :: alist2
	integer, dimension(2,3) :: a
        integer, dimension(6,4) ::X
        !Ntime+1,Nm
!subroutine calls below should be modified as required and variables should be declared above as needed.

        !call generate2(8,6,4000,qmax)
                !call adjacency_matrix()
        !(/(/1,2,3/),(/4,5,6/),(/7,8,9/)/)
        !call stats(5,3,400,100,qnetm,qmaxm,qvarm)
        !print *, 'hello'
        !call test_stats_omp(5,3,400,100,4,walltime)
       ! call adjacency_matrix()
        !call connectivity(c)
        !call vary_connectivity(8,6,4000,carray)
        a(1,1)=1
        a(1,2)=1
        a(2,1)=2
        a(2,2)=3
        a(2,3)=5
        a(1,3)=2
        X=0
        !print *, a
        !call adjacency_list([2,2],a,alist1,alist2)
       ! rwnet(Ntime,Nm,X0,N0,L,Nt,isample,X,XM)
        call rwnet(1,3,6,7,4,1,3,X,2)
        
end program test
