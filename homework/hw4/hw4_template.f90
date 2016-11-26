!Template code for hw4
!Should be compiled with network.f90
module netstats
	use network
	use omp_lib

	contains

subroutine stats(n0,l,nt,m,qnetm,qmaxm,qvarm)
	!Input: 
	!n0,nl,nt: recursive network model parameters
	!m: number of network realizations to compute
	!Output:
	!qnetm: node lists for all m networks
	!qmaxm,qvarm: ensemble averages of max(q) and var(q)
	implicit none
	integer, intent(in) :: n0,l,nt,m
	integer, dimension(n0+nt,m), intent(out) :: qnetm
	integer, intent(out) :: qmaxm
	real(kind=8), intent(out) :: qvarm


end subroutine stats


subroutine stats_omp(n0,l,nt,m,numThreads,qnetm,qmaxm,qvarm)
	!Input: 
	!n0,nl,nt: recursive network model parameters
	!m: number of network realizations to compute
	!numThreads: number of threads for parallel computation
	!Output:
	!qnetm: node lists for all m networks
	!qmaxm,qvarm: ensemble averages of max(q) and var(q)
	implicit none
	integer, intent(in) :: n0,l,nt,m,numThreads
	integer, dimension(n0+nt,m), intent(out) :: qnetm
	integer, intent(out) :: qmaxm
	real(kind=8), intent(out) :: qvarm


end subroutine stats_omp


subroutine test_stats_omp(n0,l,nt,m,numThreads,walltime)
	!Input: same as stats_omp
	!Output: walltime: time for 100 cals to stats_par
	implicit none
	integer, intent(in) :: n0,l,nt,m,numThreads
	real(kind=8), intent(out) :: walltime


end subroutine test_stats_omp


end module netstats
