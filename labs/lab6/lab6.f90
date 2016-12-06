!Use dgesv, to solve AX=B (where A is a N x N matrix of random variabels) and output timing information
program dgesv_time
    implicit none
    !declare variables
    integer :: N,NRHS,LDA,LDB,INFO
    double precision, allocatable, dimension(:,:) :: A,B,Atemp,Btemp,x
    integer, allocatable, dimension(:) :: IPIV
    real(kind=8) :: t1,t2,timeL,timeM


    !set dimensions of A,B
    open(unit=11,file='data.in')
    read(11,*) N

    NRHS = 1
    LDA = N
    LDB = N

    !allocate arrays, create matrices A and B
    allocate(A(LDA,N),B(LDB,NRHS),IPIV(N))
    allocate(Atemp(size(A,1),size(A,2)),Btemp(LDB,NRHS))

    !Generate random matrices
    call random_number(A)
    call random_number(B)


    !solve Ax = B
    Atemp = A
    Btemp = B

    call dgesv(N, NRHS, Atemp, LDA, IPIV, Btemp, LDB, INFO)

    !extract soln from Btemp
    allocate(x(N,NRHS))
    x = Btemp(1:N,:)

    print *, 'test:',maxval(abs(matmul(A,x)-B))

    print *, 'time=',timeL,timeM

end program dgesv_time


