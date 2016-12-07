!Project part 1
module rwmodule
    use network

contains



subroutine rwnet(Ntime,Nm,X0,N0,L,Nt,isample,X,XM)
    !random walks on a recursive network
    !Input variables:
    !Ntime: number of time steps
    !m: number of walks
    !X0: initial node, (node with maximum degree if X0=0)
    !N0,L,Nt: recursive network parameters
    !Output: X, m random walks on network each with Ntime steps from initial node, X0
    !XM: fraction of walks at initial node, computed every isample timesteps
	implicit none
    integer :: i, i1, i2, i3, i4, i5, i6, node, y
    integer, intent(in) :: Ntime,Nm,N0,L,Nt,X0,XM,isample
    real(kind=8), allocatable, dimension(:) :: prob
    integer, dimension(Ntime+1,Nm), intent(out) :: X
    integer :: qmax
    integer, dimension(N0+Nt) :: qnet
    integer, dimension(N0+L*Nt,2) :: enet 
    integer, allocatable, dimension(:) :: alist1
    integer, allocatable, dimension(:) :: alist2
	
    !do i=1,Nm
       call generate(N0,L,Nt,qmax,qnet,enet)
       allocate(alist1(size(enet,1)*2),alist2(size(qnet)))
       alist1=0
       alist2=0
       call adjacency_list(qnet,enet,alist1,alist2)
       allocate(prob(size(qnet)))
       prob=0
       node=X0
       !print *, prob
       y=(alist2(X0+1)-1)
         do i1=alist2(X0),y
             prob(alist1(i1))=prob(alist1(i1))+1
         end do
         !print *, prob
         
         prob=prob/sum(prob)
         print *, sum(prob)
         !print *, alist1
         !print *, size(alist1), size(prob)
    !end do
end subroutine rwnet


subroutine rwnet_omp(Ntime,Nm,X0,N0,L,Nt,isample,numthreads,X,XM)
    !parallelized version of rwnet, parallel regions should
    !use numthreads threads
	implicit none
    integer, intent(in) :: Ntime,Nm,N0,L,Nt,X0,XM,isample,numthreads
    integer, dimension(Ntime+1,Nm), intent(out) :: X



end subroutine rwnet_omp


end module rwmodule
