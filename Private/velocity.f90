subroutine velocity_update()
    USE aparallel
    USE aVariables
    implicit none
    integer :: i,j
    do i = 1,nx
        do j = 1,ny
            u(j,i) = un(j,i)    &
                -un(j,i)*dt/dx*(un(j,i)-un(j,i-1))  &
                -vn(j,i)*dt/dy*(un(j,i)-un(j-1,i))  &
                -.5d0*dt/dx/rho*(pn(j,i+1)-pn(j,i-1))   &
                +nu*dt/dx/dx*(un(j,i+1)-2*un(j,i)+un(j,i-1))    &
                +nu*dt/dy/dy*(un(j+1,i)-2*un(j,i)+un(j-1,i))
            v(j,i) = vn(j,i)    &
                -un(j,i)*dt/dx*(vn(j,i)-vn(j,i-1))  &
                -vn(j,i)*dt/dy*(vn(j,i)-vn(j-1,i))  &
                -.5d0*dt/dy/rho*(pn(j+1,i)-pn(j-1,i))   &
                +nu*dt/dx/dx*(vn(j,i+1)-2*vn(j,i)+vn(j,i-1))    &
                +nu*dt/dy/dy*(vn(j+1,i)-2*vn(j,i)+vn(j-1,i))
        end do
    end do
end subroutine
subroutine velocity_interface()
    USE aparallel
    USE aVariables
    implicit none
    un(1:ny,1:nx) = u(1:ny,1:nx)
    vn(1:ny,1:nx) = v(1:ny,1:nx)
    !   U Interface
    call interface(u,un,nx,ny)
    !   V Interface
    call interface(v,vn,nx,ny)
    call MPI_Barrier(MPI_COMM_WORLD,mpicode)
    if (my_up .le. mpirank) then
        ! Upper Boundary
        un(ny+1,:) = 1.d0
        vn(ny+1,:) = 0.d0
    end if
    if (my_down .ge. mpirank) then
        ! Lower Boundary
        un(0,:) = 0.d0
        vn(0,:) = 0.d0
    end if
    if (my_right .le. mpirank) then
        ! right Boundary
        un(:,nx+1) = 0.d0
        vn(:,nx+1) = 0.d0
    end if
    if (my_left .ge. mpirank) then
        ! left Boundary
        un(:,0) = 0.d0
        vn(:,0) = 0.d0
    end if  
end subroutine