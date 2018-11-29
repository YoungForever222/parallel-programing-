program cavity_flow
    USE aparallel
    USE aVariables
    implicit none
    call init_parallel()
    call check_parallel()
    call init_Variales()

    call end_Variables()
    call end_parallel()
end program
subroutine update()
    USE aparallel
    USE aVariables
    implicit none
end subroutine
subroutine pressure_poisson()
    USE aparallel
    USE aVariables
    implicit none
end subroutine
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
subroutine build_up_b()
    USE aparallel
    USE aVariables
    implicit none
    integer :: i,j
    do i = 1,nx
        do j = 1,ny
            b(j,i) = .5d0/dt*((un(j,i+1)-un(j,i-1))/dx + (vn(j+1,i)-vn(j-1,i))/dy)  &
                    -.25d0/dx/dx*(un(j,i+1)-un(j,i-1))**2.d0
                    -.5d0/dx/dy*(un(j+1,i)-un(j-1,i))*(vn(j,i+1)-vn(j,i-1))
                    -.25d0/dy/dy*(vn(j+1,i)-vn(j-1,i))**2.d0
        end do
    end do 
end subroutine
subroutine set_boundary()
    USE aparallel
    USE aVariables
end subroutine

    
