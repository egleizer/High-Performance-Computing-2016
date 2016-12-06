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
	real(kind=8) :: x, var, partial_var, threadID
	integer :: i1,i2, qmax
	integer, dimension(n0+nt) :: qnet
	integer, dimension(n0+l*nt,2) :: enet 
	
	qnetm(:,:)=0
	qmaxm=0
	!$OMP parallel do
	do i1=1,m
	   call generate(n0,l,nt,qmax,qnet,enet)
	   qnetm(:,i1)=qnet
	   
	   if (qmax>qmaxm) then
	        qmaxm=qmax
	   end if
	   x = sum(qnet)/size(qnet)
	   
	   var=0.d0
	   partial_var=0.d0
           !$OMP parallel firstprivate(partial_var),private(threadID)
           !$OMP do
	   do i2=1,size(qnet)
                partial_var=partial_var + (qnet(i2)-x)**2
           end do
           !$OMP end do
           !$OMP critical
           threadID = omp_get_thread_num()
           !print *, 'Thread number:',threadID, 'var=', var
           var = var + partial_var
           !$OMP end critical   
           !$OMP end parallel
           
	   qvarm=qvarm + var/size(qnet)

	end do
	!$OMP end parallel do
	qvarm=qvarm/m
        
        !print *, qvarm


end subroutine stats_omp
subroutine test_stats_omp(n0,l,nt,m,numThreads,walltime)
	!Input: same as stats_omp
	!Output: walltime: time for 100 cals to stats_par
	implicit none
	integer, intent(in) :: n0,l,nt,m,numThreads
	real(kind=8), intent(out) :: walltime
        integer, dimension(n0+nt,m) :: qnetm
	integer  :: qmaxm, i, clock_t1, clock_t2, clock_rate
	real(kind=8) :: qvarm,t
	t=0.d0
	if (numThreads==0) then 
	   do i=1,100
	     call system_clock(clock_t1)
	     call stats(n0,l,nt,m,qnetm,qmaxm,qvarm)
	     call system_clock(clock_t2,clock_rate)
	     t=t+dble(clock_t2-clock_t1)/dble(clock_rate)
	   end do
	else 
	   do i=1,100
	     call system_clock(clock_t1)
	     call stats_omp(n0,l,nt,m,numThreads,qnetm,qmaxm,qvarm)
	     call system_clock(clock_t2,clock_rate)
	     t=t+dble(clock_t2-clock_t1)/dble(clock_rate)
	   end do
	end if
	t=t/100
        print *, t
end subroutine test_stats_omp



  var=0.d0
	   partial_var=0.d0
           !$OMP parallel firstprivate(partial_var),private(threadID)
           !$OMP do
	   do i2=1,size(qnet)
                partial_var=partial_var + (qnet(i2)-x)**2
           end do
           !$OMP end do
           !$OMP critical
           threadID = omp_get_thread_num()
           !print *, 'Thread number:',threadID, 'var=', var
           var = var + partial_var
           !$OMP end critical   
           !$OMP end parallel
           