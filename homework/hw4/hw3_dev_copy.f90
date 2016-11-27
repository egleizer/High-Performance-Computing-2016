module network2
    implicit none
    ! Declaring all the dummy variables needed throughout the module
    integer :: i1, i2, i3, i4, i5, i6, i7, i8, i9, i10, i11, i12, i13, i14
    integer, allocatable, dimension(:) :: qnet
    integer, allocatable, dimension(:,:) :: anet, enet
    save
contains


subroutine generate2(N0,L,Nt,qmax)
	implicit none
	integer, intent(in) :: N0,L,Nt
	integer, intent(out) :: qmax
	real(kind=8) :: u, k
	real(kind=8), dimension(N0+Nt) :: probs
	real(kind=8), dimension(L) :: t
	!Allocating qnet, anet, enet with appropriate size
	allocate(qnet(N0+Nt),anet((N0+Nt),(N0+Nt)),enet(N0+(L*Nt),2))
	!Assigning the initial connections between N0 nodes
        enet(:,:)=0
        enet(N0,1)=1
        enet(N0,2)=dble(N0)
        qnet(:)=0.d0
     
        probs(1)=0
        !Loop for assigning the initial nodes connections
        do i1=1,N0-1
            enet(i1,1)=dble(i1)
            enet(i1,2)=dble(i1+1)
        end do
        !Loop for assigning the initial degree of nodes
         do i2=1,N0
             qnet(i2)=2
         end do
        !Loop over all new nodes        
        do i5=N0+1,(N0+Nt)
               !Loop to assign probability values to all the nodes up to the
               !one in which's loop we are
               do i6=1,(i5-1)
               probs(i6)=qnet(i6)/dble(sum(qnet))
               end do
           !Loop for creating a matrix of random numbers of size L
           do i8=1,L
            call random_number(u)
            t(i8)=u
           end do
               !Loop over L connections we need to make
               do i9=1,L
                 k=probs(1)
                 i7=1.d0    
                    !Loop to find in which part of [0,1) our cummulitative 
                    !probability lies to assign it no one of the nodes
                    do while (k<t(i9)) 
                       i7=i7+1.d0
                       k=k+probs(i7)
                    end do
                 !Assign new established connection to enet   
                 i11=N0+(i5-1-N0)*L+i9
                 enet(i11,1)=i5
                 enet(i11,2)=i7
                 !Assign new established connection to qnet
                 qnet(i7)=qnet(i7)+1
                 qnet(i5)=qnet(i5)+1
               end do
               probs(:)=0
           t(:)=0
        end do
        !Loop to find the maximum value of qnet     
        qmax=qnet(1)
        do i10=1,(N0+Nt)
           if (qnet(i10)>qmax) then
               qmax=qnet(i10)
           end if
        end do
        !print *, qnet
end subroutine generate2


end module network2
