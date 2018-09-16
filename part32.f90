!-------------------------------
module syncmodule
	implicit none
	complex(kind=16), parameter :: ii=cmplx(0.0,1.0) !ii = sqrt(-1)
        integer :: ntotal, a !total number of oscillators, 
	real(kind=8) :: c,mu,sigma !coupling coefficient, mean, std
	save
end module syncmodule
!-------------------------------

program sync_mpi
    use mpi
    use syncmodule
    implicit none
    integer :: i1,j1
    integer :: nt !number of time steps, number of oscillators
    real(kind=8) :: dt,pi !time step
    integer :: myid, numprocs, ierr
    real(kind=8), allocatable, dimension(:) :: f0,w,f ! initial condition, frequencies, solution
    real(kind=8), allocatable, dimension(:) :: order !order parameter

 ! Initialize MPI
    call MPI_INIT(ierr)
    call MPI_COMM_SIZE(MPI_COMM_WORLD, numprocs, ierr)

!gather input
    open(unit=10,file='data.in')
        read(10,*) ntotal
        read(10,*) nt
        read(10,*) dt
        read(10,*) c
        read(10,*) sigma
        read(10,*) a
    close(10)

    allocate(f0(ntotal),f(ntotal),w(ntotal),order(nt))
	

!generate initial condition
    pi = acos(-1.d0)
    call random_number(f0)
    f0 = f0*2.d0*pi


!generate frequencies
    mu = 1.d0       
    call random_normal(ntotal,w)
    w = sigma*w+mu    
    
!compute solution
    call euler_mpi(MPI_COMM_WORLD,numprocs,ntotal,0.d0,f0,w,dt,nt,f,order)


!output solution (after completion of gather in euler_mpi)
       call MPI_COMM_RANK(MPI_COMM_WORLD, myid, ierr)
       if (myid==0) then
        open(unit=11,file='theta.dat')
        do i1=1,ntotal
            write(11,*) f(i1)
        end do
        close(11)
        
        open(unit=12,file='order.dat')
        do i1=1,nt
	    write(12,*) order(i1)
	end do
	close(12)
    end if
    !can be loaded in python with: f=np.loadtxt('theta.dat')
   
    call MPI_FINALIZE(ierr)
end program sync_mpi


!subroutine output



subroutine euler_mpi(comm,numprocs,n,t0,y0,w,dt,nt,y,order)
    !explicit Euler method, parallelized with mpi
    !input: 
    !comm: MPI communicator
    !numprocs: total number of processes
    !n: number of oscillators
    !t0: initial time
    !y0: initial phases of oscillators
    !w: array of frequencies, omega_i
    !dt: time step
    !nt: number of time steps
    !output: y, final solution
    !order: order at each time step
    use mpi
    use syncmodule
    implicit none
    integer, intent (in) :: n,nt
    real(kind=8), dimension(n), intent(in) :: y0,w
    real(kind=8), intent(in) :: t0,dt
    real(kind=8), dimension(nt), intent(out) :: order
    real(kind=8), dimension(n), intent(out) :: y
    real(kind=8) :: t
    integer :: i1,k,istart,iend
    integer :: comm, myid,ierr,numprocs, sender, receiver, nn
    real(kind=8), allocatable, dimension(:) :: ylocal, Rpart, RHS
    complex(kind=8) :: e
    integer :: i2,i3,i4
    integer, allocatable, dimension(:) :: Nper_proc, disps
    integer, dimension(MPI_STATUS_SIZE) :: status
    complex(kind=16) :: summ
    complex(kind=8) :: summl, summt

    call MPI_COMM_RANK(comm, myid, ierr)
    print *, 'start euler_mpi, myid=',myid

    !set initial conditions
    y = y0
    t = t0

    !generate decomposition and allocate sub-domain variables
    call mpe_decomp1d(size(y),numprocs,myid,istart,iend)
    print *, 'istart,iend,threadID=',istart,iend,myid
    !define nn
    nn = iend - istart + 1
    !Allocate the used variables
    allocate(ylocal(nn+2*a), Rpart(nn+2*a),RHS(nn+2*a))
    
    !Loop to assign the first values used
    if (myid==numprocs-1) then
       ylocal(1:(nn+a))=y((ntotal-nn-a+1):ntotal)
       ylocal((nn+a+1):(nn+2*a))=y(1:a)  
    else if (myid==0) then
       ylocal(1:a)=y((ntotal-a+1):ntotal)
       ylocal((a+1):(nn+2*a))=y(1:(nn+a))
    else 
       ylocal=y((istart-a):(iend+a))   
    end if
       
    !time marching
    do k = 1,nt  
        !Loops to assign what is received and what is sent
        if (myid<numprocs-1) then
             receiver = myid+1
        else
             receiver = 0
        end if
 
        if (myid>0) then
             sender = myid-1
        else
             sender = numprocs-1
        end if
        
        !Sending and receiving
        call MPI_SEND(ylocal((nn+1):(nn+a)),a,MPI_DOUBLE_PRECISION,receiver,0, MPI_COMM_WORLD,ierr)
        call MPI_RECV(ylocal(1:a),a,MPI_DOUBLE_PRECISION,sender,MPI_ANY_TAG,MPI_COMM_WORLD,status,ierr)
       
        !Second pair of loops to assign what is sent and what is received
        if (myid<(numprocs-1)) then
             sender = myid+1
        else
             sender = 0
        end if
               
        if (myid>0) then
             receiver = myid-1
        else
             receiver = numprocs-1
        end if
        
        !Sending and receiving
        call MPI_SEND(ylocal((a+1):(2*a)),a,MPI_DOUBLE_PRECISION,receiver,0, MPI_COMM_WORLD,ierr)
        call MPI_RECV(ylocal((nn+a+1):(nn+2*a)),a,MPI_DOUBLE_PRECISION,sender,MPI_ANY_TAG,MPI_COMM_WORLD,status,ierr)
       
        !Calling our RHS
        call RHS_mpi(nn,t,w(istart:iend),ylocal,Rpart)
        
        
        ylocal= ylocal + dt*Rpart
	!Computing the order, using MPI_reduced as we've not gathered yet
	
	summ=0.d0
	do i4=(1+a),(nn+a)
	   summl=summl+exp(ii*ylocal(i4))
	end do
	
	call MPI_REDUCE(summl, summt, 1, MPI_DOUBLE_COMPLEX,MPI_sum,0,MPI_COMM_WORLD,ierr)
	
	if (myid==0) then
	    order(k)=(1.d0/ntotal)*abs(summt)
	end if   
    
     end do
 
 
    print *, 'before collection',myid, maxval(abs(ylocal((a+1):(nn+a))))
   !Allocating variables for MPI_gatherv
    allocate(Nper_proc(Numprocs),disps(Numprocs))
    call MPI_GATHER(nn,1,MPI_INT,Nper_proc,1,MPI_INT,0,MPI_COMM_WORLD,ierr)
      
    if (myid==0) then
        print *,'Nper_proc=',Nper_proc
        disps(1)=0
        do i1=2,Numprocs
            disps(i1) = disps(i1-1)+Nper_proc(i1-1) 
        end do
        print *, 'disps=', disps
    end if

    !collecting ylocal from each processor onto myid=0
    call MPI_GATHERV(ylocal(a+1:nn+a),nn,MPI_DOUBLE_PRECISION,y,nper_proc,&
            disps,MPI_DOUBLE_PRECISION,0,comm,ierr)

    if (myid==0) print *, 'finished',maxval(abs(y))

end subroutine euler_mpi

!-------------------------

subroutine RHS_mpi(nn,t,w,f,rhs)
    !called by euler_mpi
    !rhs = df/dt
    use syncmodule
    implicit none
    integer, intent(in) :: nn
    integer :: i1,i2
    real(kind=8), intent(in) :: t
!dimensions of variables below must be added    
    real(kind=8), dimension(nn), intent(in) :: w 
    real(kind=8), dimension(nn+2*a), intent(in) :: f
    real(kind=8), dimension(nn+2*a), intent(out) :: rhs
    real(kind=8) :: summ
 !Loop to computr RHS
    do i1=a+1,nn+a
        summ=0.d0
        summ=sum(sin(f(i1)-f((i1-a):(i1+a))))
        rhs(i1)=w(i1-a) - (c/dble(nn))*summ
    end do

end subroutine RHS_mpi


!--------------------------------------------------------------------
!  (C) 2001 by Argonne National Laboratory.
!      See COPYRIGHT in online MPE documentation.
!  This file contains a routine for producing a decomposition of a 1-d array
!  when given a number of processors.  It may be used in "direct" product
!  decomposition.  The values returned assume a "global" domain in [1:n]
!
subroutine MPE_DECOMP1D( n, numprocs, myid, s, e )
    implicit none
    integer :: n, numprocs, myid, s, e
    integer :: nlocal
    integer :: deficit

    nlocal  = n / numprocs
    s       = myid * nlocal + 1
    deficit = mod(n,numprocs)
    s       = s + min(myid,deficit)
    if (myid .lt. deficit) then
        nlocal = nlocal + 1
    endif
    e = s + nlocal - 1
    if (e .gt. n .or. myid .eq. numprocs-1) e = n

end subroutine MPE_DECOMP1D

!--------------------------------------------------------------------

subroutine random_normal(n,rn)

! Adapted from the following Fortran 77 code
!      ALGORITHM 712, COLLECTED ALGORITHMS FROM ACM.
!      THIS WORK PUBLISHED IN TRANSACTIONS ON MATHEMATICAL SOFTWARE,
!      VOL. 18, NO. 4, DECEMBER, 1992, PP. 434-435.

!  The function random_normal() returns a normally distributed pseudo-random
!  number with zero mean and unit variance.

!  The algorithm uses the ratio of uniforms method of A.J. Kinderman
!  and J.F. Monahan augmented with quadratic bounding curves.

IMPLICIT NONE
integer, intent(in) :: n
real(kind=8), intent(out) :: rn(n)
!     Local variables
integer :: i1
REAL(kind=8)     :: s = 0.449871, t = -0.386595, a = 0.19600, b = 0.25472,           &
            r1 = 0.27597, r2 = 0.27846, u, v, x, y, q

!     Generate P = (u,v) uniform in rectangle enclosing acceptance region
do i1=1,n

DO
  CALL RANDOM_NUMBER(u)
  CALL RANDOM_NUMBER(v)
  v = 1.7156d0 * (v - 0.5d0)

!     Evaluate the quadratic form
  x = u - s
  y = ABS(v) - t
  q = x**2 + y*(a*y - b*x)

!     Accept P if inside inner ellipse
  IF (q < r1) EXIT
!     Reject P if outside outer ellipse
  IF (q > r2) CYCLE
!     Reject P if outside acceptance region
  IF (v**2 < -4.d0*LOG(u)*u**2) EXIT
END DO

!     Return ratio of P's coordinates as the normal deviate
rn(i1) = v/u
end do
RETURN

!The explanation for the part 3.6


!The way we want to choose the partition depends on how well connected the network is.
!If it has a lot of disconnected parts, we want to choose the partition which would help 
!us to get as many distinct parts as possible as then we would be able to send then with MPI
!to different receivers and work with them separately. To check that we need to take a partition 
!which is more  separable and time the parallelisation with MPI. If it gives significant improvement, 
!this means we should go on with it. If it doesn’t we would rather use OMP because this means that our 
!network is not really separable. The easier way to do it would be to look at the whole structure of the network, 
!but as we  don’t have it, we can not really say if the different in partitions one and two are significant or not. 

!In fact, to be precise, if we find that setting of the MPI is too long (as it’s usability and syntactic tend to 
!be more complicated), we should rather use OMP to test the time (choosing less separated partition) and see if we get
!significant improvement or not. This depends on the way the code will look like, which is hard to describe without seeing it.
END subroutine random_normal

