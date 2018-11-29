subroutine pressure_poisson()
    USE aparallel
    USE aVariables
    implicit none
    integer :: i,j,q
    do q = 1,nit 
        call pressure_interface()
        do i = 1,nx
            do j = 1,ny
                p(j,i) = .5d0*((pn(j,i+1)+pn(j,i-1))*dy*dy+(pn(j+1,i)+pn(j-1,i))*dx*dx)/(dx*dx+dy*dy)-b(j,i)
            end do
        end do  
    end do
    call pressure_interface()   
end subroutine
subroutine pressure_interface()
    USE aparallel
    USE aVariables
    implicit none
    pn(1:ny,1:nx) = p(1:ny,1:nx)
    !   P Interface
    call interface(p,pn,nx,ny)
    if (my_up .le. mpirank) then
        ! Upper Boundary
        pn(ny+1,:) = 0.d0
    end if
    if (my_down .ge. mpirank) then
        ! Lower Boundary
        pn(0,:) = pn(1,:)
    end if
    if (my_right .le. mpirank) then
        ! right Boundary
        pn(:,nx+1) = pn(:,nx)
    end if
    if (my_left .ge. mpirank) then
        ! left Boundary
        pn(:,0) = pn(:,1)
    end if  
end subroutine
subroutine build_up_b()
    USE aparallel
    USE aVariables
    implicit none
    integer :: i,j
    do i = 1,nx
        do j = 1,ny
            b(j,i) = .5d0/dt*((un(j,i+1)-un(j,i-1))/dx + (vn(j+1,i)-vn(j-1,i))/dy)  &
                    -.25d0/dx/dx*(un(j,i+1)-un(j,i-1))**2    &
                    -.5d0/dx/dy*(un(j+1,i)-un(j-1,i))*(vn(j,i+1)-vn(j,i-1)) &
                    -.25d0/dy/dy*(vn(j+1,i)-vn(j-1,i))**2
        end do
    end do
    b(:,:) = .5d0*rho*dx*dx*dy*dy/(dx*dx+dy*dy)*b(j,i)
end subroutine
