!Evgeniia Gleizer, CID: 00948999
module netstats
	use network
	use omp_lib
        save
contains

subroutine stats(n0,l,nt,m,qnetm,qmaxm,qvarm)
	implicit none
	integer, intent(in) :: n0,l,nt, m
	integer, dimension(n0+nt,m), intent(out) :: qnetm
	integer, intent(out) :: qmaxm
	real(kind=8), intent(out) :: qvarm
	real(kind=8) :: x, var
	integer :: i1,i2, qmax
	integer, dimension(n0+nt) :: qnet
	integer, dimension(n0+l*nt,2) :: enet 
	!Assign tqnet and qmax to 0
	qnetm(:,:)=0
	qmaxm=0
	!Loop to run all further procedures m times
	do i1=1,m
           !Calling our generate function
	   call generate(n0,l,nt,qmax,qnet,enet)
	   !Putting qnet values in qnetm
	   qnetm(:,i1)=qnet
	   !Choosing max out of all qmax for different m
	   qmaxm=max(qmax,qmaxm)
	   var=0
	   !Calculating mean 
	   x = sum(qnet)/size(qnet)
	   !Loop to calculate variance 
	   do i2=1,size(qnet)
                var=var + (qnet(i2)-x)**2
           end do
	   qvarm=qvarm + var/size(qnet)

	end do
	!Diving by m to get the mean
	qvarm=qvarm/m
    
end subroutine stats

subroutine stats_omp(n0,l,nt,m,numThreads,qnetm,qmaxm,qvarm)
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
	!Parallelisation with OMP
	!$OMP parallel do private(qmax,qnet,enet,qmaxm,var,x)
	do i1=1,m
	   call generate(n0,l,nt,qmax,qnet,enet)
	   qnetm(:,i1)=qnet
	   qmaxm=max(qmax,qmaxm)
	   x = sum(qnet)/size(qnet)
	   
	   do i2=1,size(qnet)
                var=var + (qnet(i2)-x)**2
           end do
	   qvarm=qvarm + var/size(qnet)

	end do
	!End of parallelisation
	!$OMP end parallel do
	qvarm=qvarm/m

end subroutine stats_omp

subroutine test_stats_omp(n0,l,nt,m,numThreads,walltime)
	implicit none
	integer, intent(in) :: n0,l,nt,m,numThreads
	real(kind=8), intent(out) :: walltime
        integer, dimension(n0+nt,m) :: qnetm
	integer(kind=8)  :: clock_t1, clock_t2, clock_rate
	integer :: qmaxm, i
	real(kind=8) :: qvarm
	!Assigning the time variables to 0
	walltime=0.d0
	clock_t1=0.d0
	clock_t2=0.d0
	!If loop to distinguish different amount of numThreads
	if (numThreads==1) then 
	   !Loop to calcullate array of times of stats for all m
	   do i=1,100
	     call system_clock(clock_t1)
	     call stats(n0,l,nt,m,qnetm,qmaxm,qvarm)
	     call system_clock(clock_t2,clock_rate)
	     walltime=walltime+dble(clock_t2-clock_t1)/dble(clock_rate)
	   end do
	else 
	   !Loop to calcullate array of times of statsomp for all m
	   do i=1,100
	     call system_clock(clock_t1)
	     call stats_omp(n0,l,nt,m,numThreads,qnetm,qmaxm,qvarm)
	     call system_clock(clock_t2,clock_rate)
	     walltime=walltime+dble(clock_t2-clock_t1)/dble(clock_rate)
	   end do
	end if
	!Calculating the average of walltime
	walltime=walltime/100
        print *, walltime
end subroutine test_stats_omp


end module netstats
