module aparallel
    USE MPI
    implicit none
    !include 'mpif.h'
    integer,parameter :: Master = 0
    integer,parameter :: MPI_DP = MPI_DOUBLE_PRECISION
    logical :: root=.false.
    integer :: mpisize,mpirank
    integer :: mpicode
    integer :: status(MPI_STATUS_SIZE)
    character(len=10) :: nodename
    integer :: resultlen
    integer :: requestA,requestB,requestC
    integer :: my_row,my_col,my_up,my_down
end module
module Matrix
    implicit none
    integer,parameter :: r=10,s=4,t=10    ! A(r,s)*B(s,t) = C(r,t)
    real(8),allocatable,dimension(:,:) :: X,Y,Z,Serial
    real(8),allocatable,dimension(:,:) :: A,B,C
    real(8),allocatable,dimension(:,:) :: A_temp,B_temp
    integer,allocatable,dimension(:)   :: pos
    integer :: i,j,k
    integer :: p_sqrt
    integer :: dim_real(3) , dim_ext(3) , dim_sub(3)
    contains
        subroutine init()
            !   Extend to a matrix that can be divided equally
            dim_real(1) = r
            dim_real(2) = s
            dim_real(3) = t
            do i = 1,3
                j = 0
                do while(j*p_sqrt < dim_real(i)) 
                    j = j + 1
                end do
                dim_ext(i) = j*p_sqrt
                dim_sub(i) = j
            end do
            allocate(pos(1:p_sqrt))
        end subroutine
        subroutine init_random()
            implicit none
            real(8) :: temp
            call random_seed()
            allocate(X(dim_ext(1),dim_ext(2)),Y(dim_ext(2),dim_ext(3)),Z(dim_ext(1),dim_ext(3)))
            allocate(Serial(dim_ext(1),dim_ext(3)))
            X(:,:) = 0.d0
            Y(:,:) = 0.d0
            Z(:,:) = 0.d0
            do j = 1,s
                do i = 1,r
                    call random_number(temp)
                    X(i,j) = i
                end do
                do k = 1,t
                    call random_number(temp)
                    Y(j,k) = k
                end do
            end do
            open(unit=6,FILE='com.log',STATUS='REPLACE')
            write(6,*) "By MATMUL = "
            Serial = MATMUL(X,Y)
            do i = 1,r
                write(6,*) Serial(i,1:t)
            end do
        end subroutine
        !subroutine init_read()
        !    implicit none
        !end subroutine
        subroutine init_sub()
            implicit none
            allocate(A(dim_sub(1),dim_sub(2)),B(dim_sub(2),dim_sub(3)),C(dim_sub(1),dim_sub(3)))
            allocate(A_temp(dim_sub(1),dim_sub(2)),B_temp(dim_sub(2),dim_sub(3)))
            A(:,:) = 0.d0
            B(:,:) = 0.d0
            C(:,:) = 0.d0
            A_temp(:,:) = 0.d0
            B_temp(:,:) = 0.d0
        end subroutine
        subroutine output()
            open(unit=6,FILE='res.log',STATUS='REPLACE')
            write(6,*) "X = "
            do i = 1,r
                write(6,*) X(i,1:s)
            end do
            write(6,*) "Y = "
            do j = 1,s
                write(6,*) Y(j,1:t)
            end do
            write(6,*) "Z = X * Y = "
            do i = 1,r
                write(6,*) Z(i,1:t)
            end do
        end subroutine
        subroutine release()
            if(root) then
                deallocate(X)
                deallocate(Y)
                deallocate(Z)
                deallocate(Serial)
            end if
            deallocate(A)
            deallocate(B)
            deallocate(C)
            deallocate(A_temp)
            deallocate(B_temp)
            deallocate(pos)
end module
subroutine Scatter()
    Use aparallel
    Use Matrix
    implicit none
    if (root) then
        do i = 0,p_sqrt-1
            do j = 0,p_sqrt-1
                A_temp(:,:) = X(i*dim_sub(1)+1:(i+1)*dim_sub(1),j*dim_sub(2)+1:(j+1)*dim_sub(2))
                B_temp(:,:) = Y(i*dim_sub(2)+1:(i+1)*dim_sub(2),j*dim_sub(3)+1:(j+1)*dim_sub(3))
                if(i.eq.0 .and. j.eq.0) then 
                    A = A_temp
                    B = B_temp
                else
                    call MPI_ISEND(A_temp,dim_sub(1)*dim_sub(2),MPI_DP,i*p_sqrt+j,0,MPI_COMM_WORLD,requestA,mpicode)
                    call MPI_ISEND(B_temp,dim_sub(2)*dim_sub(3),MPI_DP,i*p_sqrt+j,1,MPI_COMM_WORLD,requestB,mpicode)
                end if
            end do
        end do
    else
        call MPI_IRECV(A,dim_sub(1)*dim_sub(2),MPI_DP,0,0,MPI_COMM_WORLD,requestA,mpicode)
        call MPI_IRECV(B,dim_sub(2)*dim_sub(3),MPI_DP,0,1,MPI_COMM_WORLD,requestB,mpicode)
    end if
    call MPI_WAIT(requestA,status,mpicode)
    call MPI_WAIT(requestB,status,mpicode)
    if (root) then
        write(*,*) "Scatter Over"
    end if 
    !write(*,*) my_row,my_col,"=",A,B
end subroutine


program fox
    Use aparallel
    Use Matrix
    implicit none
    !parameters
    
    !Variables
    integer :: ii
    integer :: p
    call MPI_INIT( mpicode )
    call MPI_COMM_RANK( MPI_COMM_WORLD, mpirank, mpicode )
    call MPI_COMM_SIZE( MPI_COMM_WORLD, mpisize, mpicode )
    p = mpisize
    p_sqrt = ANINT(sqrt(float(mpisize)))
    if (p_sqrt*p_sqrt .ne. p) then
        write(*,*) "mpisize is not a perfect square"
        call MPI_FINALIZE(mpicode)
        stop
    end if
    my_row = mpirank / p_sqrt
    my_col = mod(mpirank,p_sqrt)
    my_up  = mod((my_row - 1) + p_sqrt,p_sqrt)
    my_down= mod(my_row + 1,p_sqrt)
    call init()
    if (mpirank .eq. Master) then
        root = .true.
        call init_random()
    end if
    call init_sub()
    call Scatter()
    !!Start!!
    !   A_i,i
    do i = 0,p_sqrt-1
        pos(i) = i
    end do
    do ii = 0,p_sqrt-1 ! turn
        write(*,*) "ii=",ii
        !   Step 1
        if(pos(my_row) .eq. my_col) then
            A_temp = A
            do j = 0,p_sqrt-1
                if (j .eq. pos(my_row)) cycle
                call MPI_SEND(A,dim_sub(1)*dim_sub(2),MPI_DP,my_row*p_sqrt+j,0,MPI_COMM_WORLD,mpicode)
            end do
        else
            call MPI_RECV(A_temp,dim_sub(1)*dim_sub(2),MPI_DP,my_row*p_sqrt+pos(my_row),0,MPI_COMM_WORLD,status,mpicode)
        end if
        call MPI_Barrier(MPI_COMM_WORLD,mpicode)
        !write(*,*) "Step 1 over"
        !   Step 2
        C = C + MATMUL(A_temp,B)
        !write(*,*) "ii=",ii,"mpirank",mpirank,"C",C
        !write(*,*) "Step 2 over"
        !   Step 3
        call MPI_SENDRECV(B,dim_sub(2)*dim_sub(3),MPI_DP,my_up*p_sqrt+my_col,my_col,    &
                        B_temp,dim_sub(2)*dim_sub(3),MPI_DP,my_down*p_sqrt+my_col,my_col,   &
                        MPI_COMM_WORLD,status,mpicode)
        call MPI_Barrier(MPI_COMM_WORLD,mpicode)
        B = B_temp
        !write(*,*) "Step 3 over"
        !   Step 4
        do i = 0,p_sqrt-1
            pos(i) = mod(pos(i) + 1,p_sqrt)
        end do
        !write(*,*) "Step 4 over"
    end do
    !! End !!
    !! Gather and Output !!
    if (root) then
        do i = 0,p_sqrt-1
            do j = 0,p_sqrt-1
                if(i.eq.0 .and. j.eq.0) then 
                    Z(i*dim_sub(1)+1:(i+1)*dim_sub(1),j*dim_sub(3)+1:(j+1)*dim_sub(3)) = C
                else
                    call MPI_RECV(Z(i*dim_sub(1)+1:(i+1)*dim_sub(1),j*dim_sub(3)+1:(j+1)*dim_sub(3)),  &
                    dim_sub(1)*dim_sub(3),MPI_DP,i*p_sqrt+j,i*p_sqrt+j,MPI_COMM_WORLD,status,mpicode)
                end if
            end do
        end do
    else
        call MPI_SEND(C,dim_sub(1)*dim_sub(3),MPI_DP,0,mpirank,MPI_COMM_WORLD,mpicode)
    end if
    call MPI_Barrier(MPI_COMM_WORLD,mpicode)
    if(root) call output()
    call MPI_FINALIZE(mpicode)
end program
        