module aparallel
    USE MPI
    !include 'mpif.h'
    integer,parameter :: Master = 0
    integer,parameter :: MPI_DP = MPI_DOUBLE_PRECISION
    logical :: root=.false.
    integer :: mpisize,mpirank,p_sqrt
    integer :: mpicode
    integer :: status(MPI_STATUS_SIZE)
    character(len=10) :: nodename
    integer :: resultlen
    integer :: requestu,requestv,requestp
    integer :: my_row,my_col,my_up,my_down,my_left,my_right
contains
    subroutine init_parallel()
        implicit none
        call MPI_INIT( mpicode )
        call MPI_COMM_RANK( MPI_COMM_WORLD, mpirank, mpicode )
        call MPI_COMM_SIZE( MPI_COMM_WORLD, mpisize, mpicode )
        if(mpirank .eq. Master) root = .true.
    end subroutine
    subroutine end_parallel()
        call MPI_FINALIZE(mpicode)
        stop
    end subroutine
end module
