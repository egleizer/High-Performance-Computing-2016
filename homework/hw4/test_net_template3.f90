program test_network2
	use netstats
	implicit none
	integer, dimension(405,100) :: qnetm
	integer :: qmaxm, qmax, numThreads
	real(kind=8) :: qvarm, walltime
	
!subroutine calls below should be modified as required and variables should be declared above as needed.

        !call generate2(8,6,4000,qmax)
                !call adjacency_matrix()
                !call output_anet()
                !call output_enet()
                !call output_qnet()
                !call connectivity()
        
        !call stats(5,3,400,100,qnetm,qmaxm,qvarm)
        print *, 'hello'
        call test_stats_omp(5,3,400,100,4,walltime)
       ! call adjacency_matrix()
        !call connectivity(c)
        !call vary_connectivity(8,6,4000,carray)
        
end program test_network2