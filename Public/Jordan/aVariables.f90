module aVariables
    USE aparallel
    !Parameters
    implicit none
    integer,parameter :: n_real=4
    integer:: m,n
    real(8),parameter :: epsilon = 1d-6
    real(8) :: max = epsilon + 1.d0
    real(8),allocatable,dimension(:,:) :: A
    real(8),allocatable,dimension(:) :: b
    real(8),allocatable,dimension(:) :: x,x_1,b_sub
    real(8),allocatable,dimension(:,:) :: A_sub
contains
    subroutine init_Variables()
        implicit none
        integer i,j
        real(8) :: temp
        j = 0
        do while(j*mpisize < n_real)
            j = j + 1
        end do
        m = j
        n = j*mpisize
        if (root) then
            !   Extend to a matrix that can be divided equally
            allocate(A(0:n-1,0:n-1),b(0:n-1))
            !   Initialize the Matrix
            !   diagonal matrix
            A(:,:) = 0.d0
            do i=0,n-1
                A(i,i) = 1.d0
                b(i) = 1.d0
            end do
            if(n_real.eq.4) then
                !   现代数值计算方法 P56-例1
                A(0,0) = 5.d0
                A(1,0) = 2.d0
                A(2,0) = 1.d0
                A(3,0) = -1.d0

                A(0,1) = 1.d0
                A(1,1) = 8.d0
                A(2,1) = -2.d0
                A(3,1) = 3.d0

                A(0,2) = -1.d0
                A(1,2) = 1.d0
                A(2,2) = -4.d0
                A(3,2) = 2.d0

                A(0,3) = -2.d0
                A(1,3) = 3.d0
                A(2,3) = -1.d0
                A(3,3) = 7.d0
                
                b(0) = -2.d0
                b(1) = -6.d0
                b(2) = 6.d0
                b(3) = 12.d0
            else
                call random_seed()
                do i = 0,n_real-1
                    do j = 0,n_real-1
                        call random_number(temp)
                        A(i,j) = temp
                    end do
                    call random_number(temp)
                    b(i) = temp
                end do
            end if
        end if
        !   For all processes
        allocate(x(0:n-1),x_1(0:n-1))
        x(:) = 0.d0
        x_1(:) = 0.d0
    end subroutine
    subroutine Scatter()
        implicit none
        integer :: i,j
        real(8) :: temp_sub(0:n-1,0:m-1),temp_A(0:n-1,0:n-1),temp_b(0:n-1)
        allocate(A_sub(0:m-1,0:n-1),b_sub(0:m-1))
        if(root) then
            do i = 0,mpisize-1
                do j = 0,m-1
                    temp_A(:,i*m+j) = A(i+j*mpisize,:)
                    temp_b(i*m+j) = b(i+j*mpisize)
                end do
            end do
        end if
        call MPI_Scatter(temp_A,n*m,MPI_DP,temp_sub,n*m,MPI_DP,0,MPI_COMM_WORLD,mpicode)
        call MPI_Scatter(temp_b,m  ,MPI_DP,b_sub   ,m  ,MPI_DP,0,MPI_COMM_WORLD,mpicode)
        A_sub = TRANSPOSE(temp_sub)
        write(*,*) "mpirank=",mpirank,"A_sub",A_sub,"b=",b_sub
        call MPI_Barrier(MPI_COMM_WORLD,mpicode)
    end subroutine
    subroutine Output()
        implicit none
        integer:: i,j
        if(root) then
            open(unit=6,file="./res.dat")
            write(6,*) "Ax = b"
            write(6,*) "A"
            do i = 0,n_real-1
                write(6,*) A(i,:)
            end do
            write(6,*) "b"
            do i = 0,n_real-1
                write(6,*) b(i)
            end do
            write(6,*) "x"
            do i = 0,n_real-1
                write(6,*) x(i)
            end do
            close(6)
        end if
        write(*,*) "check the res.dat for the result"
    end subroutine
end module

