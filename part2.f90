module flunet
        !use omp_lib
        !add variables as needed
	!save
contains


subroutine rhs(n,y,t,a,b0,b1,g,k,w,P,dy)
    implicit none
    !Return RHS of network flu model
    !input: 
    !n: total number of nodes
    !y: S,E,C
    !t: time
    !a,b0,b1,g,k,w: model parameters
    !output: dy, RHS
    integer, intent(in) :: n
    real(kind=8), dimension(n*3),intent(in) :: y
    real(kind=8), intent(in) :: t,a,b0,b1,g,k,w
    real(kind=8), dimension(n*3), intent(out) :: dy
    real(kind=8), dimension(n,n) :: P
    real(kind=8), dimension(n) :: S,E,C
    real(kind=8), dimension(n) :: dS, dE, dC
    integer :: i
    real(kind=8) :: pi, b
    
    
    !Assigninhg pi
    pi=4.d0*atan(1.d0)
    b = b0 + b1*(1.0+cos(2.0*pi*t))
    
    !Assigning the initial values
    S(:)=y(1:n)
    E(:)=y(n+1:2*n)
    C(:)=y(2*n+1:3*n)  
    dS(:)=0.d0
    dE(:)=0.d0
    dC(:)=0.d0

    !Creating RHS
    do i=1,n
        dS(i)= k*(1-S(i)) - b*C(i)*S(i) + w*dot_product(P(i,:),S)-w*S(i)
        dE(i)= b*C(i)*S(i) - (k+a)*E(i) + w*dot_product(P(i,:),E)-w*E(i)
        dC(i)= a*E(i) - (g+k)*C(i) + w*dot_product(P(i,:),C)-w*C(i)
    end do
   
    !Assigning output to one array
    dy(1:n)=dS(:)
    dy(n+1:2*n)=dE(:)
    dy(2*n+1:3*n)=dC(:)
        

end subroutine rhs

subroutine rhs_omp(n,y,t,a,b0,b1,g,k,w,P,numThreads,dy)
    implicit none
    !Return RHS of network flu model
    !input: 
    !n: total number of nodes
    !y: S,E,C
    !t: time
    !a,b0,b1,g,k,w: model parameters
    !output: dy, RHS
    integer, intent(in) :: n, numThreads
    real(kind=8), dimension(n*3),intent(in) :: y
    real(kind=8), intent(in) :: t,a,b0,b1,g,k,w
    real(kind=8), dimension(n*3), intent(out) :: dy
    real(kind=8), dimension(n,n) :: P
    real(kind=8), dimension(n) :: S,E,C
    real(kind=8), dimension(n) :: dS, dE, dC
    integer :: i
    real(kind=8) :: pi, b
    
    !Start parallelisation by calling omp_num_threads
    call omp_set_num_threads(numThreads)
    
    !Assigninhg pi
    pi=4.d0*atan(1.d0)
    b = b0 + b1*(1.0+cos(2.0*pi*t))
    
    !Assigning the initial values
    S(:)=y(1:n)
    E(:)=y(n+1:2*n)
    C(:)=y(2*n+1:3*n)  
    
    !Start parallelisation
    !$OMP parallel do 
    !Creating RHS
    do i=1,n
        dS(i)= k*(1-S(i)) - b*C(i)*S(i) + w*dot_product(P(i,:),S)-w*S(i)
        dE(i)= b*C(i)*S(i) - (k+a)*E(i) + w*dot_product(P(i,:),E)-w*E(i)
        dC(i)= a*E(i) - (g+k)*C(i) + w*dot_product(P(i,:),C)-w*C(i)
    end do
    !End parallelisation
    !$OMP end parallel do 
    
    !Assigning output to one array
    dy(1:n)=dS(:)
    dy(n+1:2*n)=dE(:)
    dy(2*n+1:3*n)=dC(:)

end subroutine rhs_omp

end module flunet
