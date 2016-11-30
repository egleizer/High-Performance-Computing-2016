! This code illustrates differentation of data distributed across subdomains
! Code outline:
! 1. Initialize MPI
! 2. Read Ntotal from data.in
! 3. Construct domain decomposition
! 4. Make grid and field, f=sin(2*pi*x), in the local subdomain
! 5. Compute derivative
! 6. Output error
! 7. Gather data on one processor, output it in Matlab-friendly format
!To compile this code: mpif90 -freal-4-real-8 -o gather.exe gather.f90
!To run: mpiexec -n 2 gather.exe
!----------------------------------------------------------------------
program gradient
    use mpi
	implicit none
	integer :: i1,j1,Ntotal,Nlocal,istart,iend
    double precision :: dx,buffer
    double precision, allocatable, dimension(:) :: x,f,df,exact
    integer :: myid, numprocs, ierr, sender, receiver
    integer, dimension(MPI_STATUS_SIZE) :: status
    double precision, parameter :: pi = 3.14159265358979323846

    !for gather:
    integer, allocatable, dimension(:) :: Nper_proc, disps
    double precision, allocatable, dimension(:) :: f_total,df_total,error_total


    ! Initialize MPI
    call MPI_INIT(ierr)
    call MPI_COMM_RANK(MPI_COMM_WORLD, myid, ierr)
    call MPI_COMM_SIZE(MPI_COMM_WORLD, numprocs, ierr)


	!read data from data.in
	open(unit=10, file='data.in')
        read(10,*) Ntotal
	close(10)

    !construct decompositio: nassign Nlocal points from istart to iend  to each processor
    call MPE_DECOMP1D( Ntotal, numprocs, myid, istart, iend)
    Nlocal = iend - istart + 1
    allocate(x(Nlocal),f(Nlocal+2),df(Nlocal),exact(Nlocal)) !set dimensions of x,f,df

    !make grid and field
    call make_grid(Ntotal,Nlocal,istart,iend,x)
    dx = x(2)-x(1)
    print *, 'proc', myid, ' has been assigned the interval x=', x(1),x(Nlocal)

    call make_field(Nlocal,x,f(2:Nlocal+1)) !note: f(1) and f(Nlocal+2) must be obtained from neighboring processors

!Now, compute derivative, df/dx
!-----------------------------------------------------------
!Send data at top boundary up to next processor
!i.e. send f(nlocal+1) to myid+1 and store it there as f(1)
!data from myid=numprocs-1 is sent to myid=0
!-----------------------------------------------------------
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

    call MPI_SEND(f(Nlocal+1),1,MPI_DOUBLE_PRECISION,receiver,0,MPI_COMM_WORLD,ierr)
    call MPI_RECV(f(1),1,MPI_DOUBLE_PRECISION,sender,MPI_ANY_TAG,MPI_COMM_WORLD,status,ierr)

    call MPI_BARRIER(MPI_COMM_WORLD,ierr)

!-----------------------------------------------------------
!Send data at bottom boundary down to previous processor
!i.e. send f(2) to myid-1 and store it there as f(nlocal+2)
!data from myid=0 is sent to myid=numprocs-1
!-----------------------------------------------------------
    if (myid>0) then
        receiver = myid-1
    else
        receiver = numprocs-1
    end if

    if (myid<numprocs-1) then
        sender = myid+1
    else
        sender = 0
    end if
    call MPI_SEND(f(2),1,MPI_DOUBLE_PRECISION,receiver,0,MPI_COMM_WORLD,ierr)
    call MPI_RECV(f(Nlocal+2),1,MPI_DOUBLE_PRECISION,sender,MPI_ANY_TAG,MPI_COMM_WORLD,status,ierr)
!----finished sending/receiving
    df = (f(3:Nlocal+2) - f(1:Nlocal))/(2.0*dx) !approximate derivative

    exact = 2.0*pi*cos(2.0*pi*x)

    print *, 'myid=', myid, 'error=', maxval(abs(exact-df))

!---------------------------------------------------------------------------------------

!Gather fields on myid=0, and output in Matlab-friendly format

    allocate(df_total(Ntotal),Nper_proc(Numprocs),disps(Numprocs))

    !gather Nlocal from each proc to array Nper_proc on myid=0
    call MPI_GATHER(Nlocal,1,MPI_INT,Nper_proc,1,MPI_INT,0,MPI_COMM_WORLD,ierr)
    if (myid==0) then
        print *,'Nper_proc=',Nper_proc
        disps(1)=0
        do i1=2,Numprocs
            disps(i1) = disps(i1-1)+Nper_proc(i1-1) !needed for gatherv below
        end do
        print *, 'disps=', disps
    end if

    !collect df from each processor onto myid=0
    call MPI_GATHERV(df,Nlocal,MPI_DOUBLE_PRECISION,df_total,Nper_proc, &
                disps,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)


    !output
    if (myid==0) then
        open(unit=23,file='field.dat')
        do i1=1,Ntotal
            write(23,*) dx*(i1-1),df_total(i1)
        end do
        close(23)
    end if
    !to plot in matlab: 
    !load field.dat
    !figure
    !plot(field(:,1),field(:,2))


    call MPI_FINALIZE(ierr)

end program gradient
!--------------------------------------------------------------------
!--------------------------------------------------------------------
!--------------------------------------------------------------------
!--------------------------------------------------------------------
!---------------------------------------------------------------------
!subroutine make_grid
!
!Let xtotal be a uniformly spaced grid from 0 to 1 with Ntotal+1 points
!This subroutine generates x = xtotal(istart:iend)
!---------------------------------------------------------------------
subroutine make_grid(Ntotal,Nlocal,istart,iend,x)
    implicit none
    integer :: i1
    integer, intent(in) :: Ntotal,Nlocal,istart,iend
    double precision, dimension(Nlocal), intent(out) :: x

    do i1 = istart,iend
        x(i1-istart+1) = dble(i1-1)/dble(Ntotal)
    end do

end subroutine make_grid
!--------------------------------------------------------------------

!---------------------------------------------------------------------
!subroutine make_field
!generates a simple sinusoidal function, f(X)
!input: N, size of field
!       x, grid points
!---------------------------------------------------------------------
subroutine make_field(N,x,f)
    implicit none
    integer, intent(in) :: N
    double precision, dimension(N), intent(in) :: x
    double precision, dimension(N), intent(out) :: f
    double precision, parameter :: pi = 3.14159265358979323846

    f = sin(2.0*pi*x)

end subroutine make_field

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





