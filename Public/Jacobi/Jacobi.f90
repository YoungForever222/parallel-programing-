Program Jacobi
    Use aparallel
    Use aVariables
    implicit none
    call init_parallel()
    call init_Variables()
    if (root)write(*,*) "init_Variables() over"
    call Scatter()
    if (root)write(*,*) "Scatter() over"
    call Jacobi_Method_MPI()
    if (root)write(*,*) "Jacobi_Method_MPI() over"
    call Output()
    call end_parallel()
end
subroutine Jacobi_Method_MPI()
    Use aparallel
    Use aVariables
    implicit none
    integer :: i,j,log
    real(8) :: sum,localmax,x_1_temp(m)
    log = 0
    do while(max .gt. epsilon)
        if (root) write(*,*) "log=",log,"max=",max
        log = log + 1
        do i = 1,m !(1)
            sum = 0.d0 !(1.1)
            do j = 1,n !(1.2)
                if(j .ne. mpirank*m+i) then
                    sum = sum + A_sub(i,j) * x(j)
                end if
            end do
            x_1(mpirank*m+i) = (b_sub(i)-sum) / A_sub(i,mpirank*m+i)
        end do
        localmax = abs(x_1(mpirank*m+1)-x(mpirank*m+1))
        do i = 2,m
            if(abs(x_1(mpirank*m+i)-x(mpirank*m+i)) .gt. localmax) then
                localmax = abs(x_1(mpirank*m+i)-x(mpirank*m+i))
            end if
        end do
        x_1_temp(:) = x_1(mpirank*m+1:mpirank*m+m)
        call MPI_Allgather(x_1_temp,m,MPI_DP, &
                            x,m,MPI_DP,MPI_COMM_WORLD,mpicode)
        call MPI_Allreduce(localmax,max,1,MPI_DP,MPI_MAX,MPI_COMM_WORLD,mpicode)
    end do
end subroutine