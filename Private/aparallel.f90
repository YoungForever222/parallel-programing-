module aparallel
    USE MPI
    !include 'mpif.h'
    integer,parameter :: Master = 0
    integer,parameter :: MPI_DP = MPI_DOUBLE_PRECISION
    logical :: root=.false.
    integer :: mpisize,mpirank,p,p_sqrt
    integer :: mpicode
    integer :: status(MPI_STATUS_SIZE)
    character(len=10) :: nodename
    integer :: resultlen
    integer :: my_row,my_col,my_up,my_down,my_left,my_right
contains
    subroutine init_parallel()
        call MPI_INIT( mpicode )
        call MPI_COMM_RANK( MPI_COMM_WORLD, mpirank, mpicode )
        call MPI_COMM_SIZE( MPI_COMM_WORLD, mpisize, mpicode )
        if(mpirank .eq. Master) root = .true.
    end subroutine
    subroutine check_parallel()
        p = mpisize
        p_sqrt = ANINT(sqrt(float(mpisize)))
        if (p_sqrt*p_sqrt .ne. p) then
            write(*,*) "mpisize is not a perfect square"
            call end_parallel()
            stop
        end if
        my_row  = mpirank / p_sqrt
        my_col  = mod(mpirank,p_sqrt)
        my_up   = mod((my_row - 1) + p_sqrt,p_sqrt)*p_sqrt+my_col
        my_down = mod(my_row + 1,p_sqrt)*p_sqrt+my_col
        my_left = my_row*p_sqrt + mod((my_col - 1) + p_sqrt,p_sqrt)
        my_right= my_row*p_sqrt + mod(my_col + 1,p_sqrt)
    end subroutine
    subroutine end_parallel()
        call MPI_FINALIZE(mpicode)
        stop
    end subroutine
end module
