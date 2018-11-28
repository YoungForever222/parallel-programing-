module aparallel
    USE MPI
    !include 'mpif.h'
    integer,parameter :: Master = 0
    integer,parameter :: MPI_DP = MPI_DOUBLE_PRECISION
    logical :: root=.false.
    integer :: mpisize,mpirank
    integer :: mpicode
    integer :: status(MPI_STATUS_SIZE)
    character(len=10) :: nodename
    integer :: resultlen
    end module
!   将属于同一节点的线程归为一个组
Program Main
    Use aparallel
    implicit none
    !parameters
    integer:: Use_Se_Re = 1
    !Variables
    logical:: first = .true.
    character(len=10),allocatable:: nameArray(:)
    integer:: MyWorld,SplitWorld
    integer:: group_size,Color,Key,Count,i
    integer:: request
    call MPI_INIT( mpicode )
    call MPI_Comm_dup(MPI_COMM_WORLD,MyWorld,mpicode)
    call MPI_COMM_RANK( MyWorld, mpirank, mpicode )
    call MPI_COMM_SIZE( MyWorld, mpisize, mpicode )
    call MPI_GET_PROCESSOR_NAME( nodename, resultlen, mpicode)
    allocate(nameArray(0:mpisize-1))
    if (mpirank == Master) root=.true.
    if(Use_Se_Re .eq. 1) then
        !MPI_SEND & MPI_RECV
        if (root) then
            nameArray(0) = nodename
            do i = 1,mpisize-1
                call MPI_RECV(nameArray(i),10,MPI_CHARACTER,i,i,MyWorld,status,mpicode)    
            end do
        else
            call MPI_SEND(nodename,10,MPI_CHARACTER,0,mpirank,MyWorld,mpicode)
        end if
        call MPI_BARRIER(MyWorld,mpicode)
        if (root) then
            do i = 1,mpisize-1
                call MPI_SEND(nameArray,10*mpisize,MPI_CHARACTER,i,i,MyWorld,mpicode)
            end do
        else
            call MPI_RECV(nameArray,10*mpisize,MPI_CHARACTER,0,mpirank,MyWorld,status,mpicode)
        end if
        call MPI_BARRIER(MyWorld,mpicode)
    else if (Use_Se_Re .eq. 2) then
        !MPI_ISEND & MPI_IRECV
        do i=0,mpisize-1
            call MPI_ISEND(nodename,10,MPI_CHARACTER,i,mpirank,MyWorld,request,mpicode)
            call MPI_IRECV(nameArray(i),10,MPI_CHARACTER,i,i,MyWorld,request,mpicode)
        end do
        call MPI_WAIT(request,status,mpicode)
    else
        call MPI_Allgather(nodename,10,MPI_CHARACTER,nameArray,10,MPI_CHARACTER,MyWorld,mpicode)
    end if
    !write(*,*) "mpirank=",mpirank,"nameArray=",nameArray
    call MPI_BARRIER(MyWorld,mpicode)
    Count = 0
    do i=0,mpirank
        if(nameArray(i) == nodename) then
            Count = Count+1
            if(first) then
                Color = i
                first = .false.
            end if
        end if
    end do
    Key = Count -1
    call MPI_Comm_split(MyWorld,Color,Key,SplitWorld,mpicode)
    write(*,*) "Color=",Color,"Key=",Key,"nodename=",nodename
    call MPI_FINALIZE(mpicode)
    end program


        
