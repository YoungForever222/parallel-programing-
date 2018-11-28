module aVariables
    USE aparallel
    !Parameters
    integer,parameter :: xDim = 1024 , yDim = 1024
    integer,parameter :: nt = 1000
    integer,parameter :: nit= 50
    integer,parameter :: c  = 1
    real(8),parameter :: Lref = 2.d0
    real(8),parameter :: dx = Lref / (xDim - 1)
    real(8),parameter :: dy = Lref / (yDim - 1)
    real(8),parameter :: rho = 1.d0, nu = 1.d0, dt = 1e-3
    !Variables
    integer :: nx , ny
    real(8),dimension(yDim,xDim)      :: u_all,v_all,p_all
    real(8),allocatable,dimension(:,:):: u,v,p,b
    real(8),allocatable,dimension(:,:):: un,vn,pn
contains
    subroutine init_Variables()
        
        allocate(un(0:ny+1,0:nx+1))
        allocate(vn(0:ny+1,0:nx+1))
        allocate(pn(0:ny+1,0:nx+1))
        allocate(b(ny,nx))
        allocate(u(ny,nx))
        allocate(v(ny,nx))
        allocate(p(ny,nx))
    end subroutine
    subroutine end_Variables()
        deallocate(u)
        deallocate(v)
        deallocate(p)
        deallocate(b)
        deallocate(un)
        deallocate(vn)
        deallocate(pn)
    end subroutine
    
