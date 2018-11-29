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
    real(8),allocatable,dimension(:,:):: u_all,v_all,p_all
    real(8),allocatable,dimension(:,:):: u,v,p,b
    real(8),allocatable,dimension(:,:):: un,vn,pn
contains
    subroutine init_Variables()
        nx = xDim/p_sqrt
        ny = yDim/p_sqrt
        if(nx*p_sqrt .ne. xDim .or. ny*p_sqrt .ne. yDim) then
            write(*,*) "The Fluid Domain Cannot Be Divided into same size"
            call end_parallel()
            stop
        end if
        if(root) then
            allocate(u_all(yDim,xDim))
            allocate(v_all(yDim,xDim))
            allocate(p_all(yDim,xDim))
        end if
        allocate(un(0:ny+1,0:nx+1))
        allocate(vn(0:ny+1,0:nx+1))
        allocate(pn(0:ny+1,0:nx+1))
        allocate(b(ny,nx))
        allocate(u(ny,nx))
        allocate(v(ny,nx))
        allocate(p(ny,nx))
    end subroutine
    subroutine end_Variables()
        if(root) then
            deallocate(u_all)
            deallocate(v_all)
            deallocate(p_all)
        end if
        deallocate(un)
        deallocate(vn)
        deallocate(pn)
        deallocate(u)
        deallocate(v)
        deallocate(p)
        deallocate(b)
    end subroutine
    
