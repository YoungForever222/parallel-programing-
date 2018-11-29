subroutine interface(a,an,nx,ny)
    USE aparallel
    implicit none
    integer :: nx,ny
    real(8) :: a(ny,nx),an(0:ny+1,0:nx+1)
    call MPI_Barrier(MPI_COMM_WORLD,mpicode)
    call MPI_SENDRECV(a(ny,:),nx,MPI_DP,my_up   ,mpirank,   &
                    an(0,1:nx)   ,nx,MPI_DP,my_down ,my_down ,MPI_COMM_WORLD,status,mpicode)
    call MPI_Barrier(MPI_COMM_WORLD,mpicode)
    call MPI_SENDRECV(a(1 ,:),nx,MPI_DP,my_down ,mpirank,   &
                    an(ny+1,1:nx),nx,MPI_DP,my_up   ,my_up   ,MPI_COMM_WORLD,status,mpicode)
    call MPI_Barrier(MPI_COMM_WORLD,mpicode)
    call MPI_SENDRECV(a(:,nx),ny,MPI_DP,my_right,mpirank,   &
                    an(1:ny,0)   ,ny,MPI_DP,my_left ,my_left ,MPI_COMM_WORLD,status,mpicode)
    call MPI_Barrier(MPI_COMM_WORLD,mpicode)
    call MPI_SENDRECV(a(:,1 ),ny,MPI_DP,my_left ,mpirank,   &
                    an(1:ny,nx+1),ny,MPI_DP,my_right,my_right,MPI_COMM_WORLD,status,mpicode)
    call MPI_Barrier(MPI_COMM_WORLD,mpicode)
end subroutine
