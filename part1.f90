!Project part 1
module rwmodule
    use network

contains



subroutine rwnet(Ntime,Nm,X0,N0,L,Nt,isample,X)
    !random walks on a recursive network
    !Input variables:
    !Ntime: number of time steps
    !m: number of walks
    !X0: initial node, (node with maximum degree if X0=0)
    !N0,L,Nt: recursive network parameters
    !Output: X, m random walks on network each with Ntime steps from initial node, X0
    !XM: fraction of walks at initial node, computed every isample timesteps
	implicit none
    integer :: i, i1, i2, i3, i4, i5, i6, i7, node, y, y1, y2, k,q,n,z, o
    integer, intent(in) :: Ntime,Nm,N0,L,Nt,X0,isample
    real(kind=8) :: u, cp
    real(kind=8), allocatable, dimension(:) :: prob
    integer, dimension(Ntime+1,Nm), intent(out) :: X
    integer :: qmax
    integer, dimension(N0+Nt) :: qnet, a
    integer, dimension(1) :: maxl
    integer, dimension(N0+L*Nt,2) :: enet 
    integer, dimension(size(enet,1)*2) :: alist1
    integer, dimension(size(qnet)) :: alist2
    integer, dimension(N0+Nt,N0+Nt) :: anet
    real(kind=8), dimension((Ntime+1)/isample) :: XM
    !Calling all the functions we will later use	
       call generate(N0,L,Nt,qmax,qnet,enet)
       call adjacency_list(qnet,enet,alist1,alist2)
       call adjacency_matrix(N0+Nt,enet,anet)
       allocate(prob(size(qnet)))
    
    !Loop for all the walks   
    do i=1,Nm
       X(1,i)=X0
       
       q=0
       !Loop to assign X0 to the maximum node in case it is given as 0
       if (X0==0) then 
          maxl=maxloc(qnet)
          X(1,i)=maxl(1)
       end if 

       !Loop over all steps
       do i3=2,(Ntime+1)
                    prob=0
                    !Finding the node we going to work with
                    k=X((i3-1),i)
                    !Assigning probability(to array, containing zeros)
                    prob=dble(anet(k,:))/qnet(X((i3-1),i))
                        
	            u=0
	            !Calling the random number to find probability
	            call random_number(u)
	            cp=0
                    !Loop to go over cumulitative probability and find the node for the next step
	            do i4=1,size(prob)
	               cp=cp+prob(i4)
	                if (u<=cp) exit
                    end do 
                !Assigning found node to the next element in the matrix    
                X(i3,i)=i4                   	 
	end do 
    end do 
        deallocate(prob)
        

   
     XM=0
     !Loop over all Ntimes to find Xm
     do i6=1,(Ntime+1),isample
         !Second part of the loop - over all Nm's
         do i7=1,Nm
             !Condition for adding the node to the ones we are interested in 
             if (X(i6,i7)==X0) then
              XM(i6)=XM(i6)+1  
             end if
         end do
     end do 
     XM=XM/Nm
     
    print *, X

end subroutine rwnet


subroutine rwnet_omp(Ntime,Nm,X0,N0,L,Nt,isample,numthreads,X)
    !parallelized version of rwnet, parallel regions should
    !use numthreads threads
    	implicit none
    integer :: i, i1, i2, i3, i4, i5, i6, i7, node, y, y1, y2, k,q,n,z, o
    integer, intent(in) :: Ntime,Nm,N0,L,Nt,X0,isample, numthreads
    real(kind=8) :: u, cp
    real(kind=8), allocatable, dimension(:) :: prob
    integer, dimension(Ntime+1,Nm), intent(out) :: X
    integer :: qmax
    integer, dimension(N0+Nt) :: qnet, a
    integer, dimension(1) :: maxl
    integer, dimension(N0+L*Nt,2) :: enet 
    integer, dimension(size(enet,1)*2) :: alist1
    integer, dimension(size(qnet)) :: alist2
    integer, dimension(N0+Nt,N0+Nt) :: anet
    real(kind=8), dimension((Ntime+1)/isample) :: XM
	
    call generate(N0,L,Nt,qmax,qnet,enet)
    call adjacency_list(qnet,enet,alist1,alist2)
    call adjacency_matrix(N0+Nt,enet,anet)
    allocate(prob(size(qnet)))
    ! Call the omp with the number of threads
    !$ call omp_set_num_threads(numThreads)
    
    !Parallising the loop
    !$OMP parallel do private(i,i3,prob,k,u,cp)
    !Loop for all the walks   
    do i=1,Nm
       X(1,i)=X0
       
       q=0
       !Loop to assign X0 to the maximum node in case it is given as 0
       if (X0==0) then 
          maxl=maxloc(qnet)
          X(1,i)=maxl(1)
       end if 

       !Loop over all steps
       do i3=2,(Ntime+1)
                    prob=0
                    !Finding the node we going to work with
                    k=X((i3-1),i)
                    !Assigning probability(to array, containing zeros)
                    prob=dble(anet(k,:))/qnet(X((i3-1),i))
                        
	            u=0
	            !Calling the random number to find probability
	            call random_number(u)
	            cp=0
                    !Loop to go over cumulitative probability and find the node for the next step
	            do i4=1,size(prob)
	               cp=cp+prob(i4)
	                if (u<=cp) exit
                    end do 
                !Assigning found node to the next element in the matrix    
                X(i3,i)=i4                   	 
	end do 
    end do 
    
        !Finishing parallelization
        !$OMP end parallel do 
        deallocate(prob)
   

 !Ntime+1,Nm  
     XM=0
     
     !Parallelizing the Xm loop
    !$OMP parallel do private(i6,i7)
     do i6=1,(Ntime+1),isample
         do i7=1,Nm
             if (X(i6,i7)==X0) then
              XM(i6)=XM(i6)+1  
             end if
         end do
     end do
     !Finishing parallelization
     !$OMP end parallel do  
     XM=XM/Nm
     


end subroutine rwnet_omp


end module rwmodule
