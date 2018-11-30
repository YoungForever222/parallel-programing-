module aVariables
    USE aparallel
    !Parameters
    implicit none
    integer,parameter :: xDim = 200 , yDim = 200
    integer,parameter :: nt = 1000000
    integer,parameter :: nit= 50
    integer,parameter :: ntout = 100000
    integer,parameter :: c  = 1
    real(8),parameter :: Lref = 2.d0
    real(8),parameter :: dx = Lref / (xDim + 1)
    real(8),parameter :: dy = Lref / (yDim + 1)
    real(8),parameter :: rho = 1.d0, nu = .1d0, dt = 1e-5
    !Variables
    integer :: nx , ny
    real(8),allocatable,dimension(:,:):: u_all,v_all,p_all,vor
    real(8),allocatable,dimension(:,:):: u,v,p,b
    real(8),allocatable,dimension(:,:):: un,vn,pn
    integer :: n = 0
contains
    subroutine init_Variables()
        implicit none
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
            allocate(vor(yDim,xDim))
            u_all(:,:) = 0.d0
            v_all(:,:) = 0.d0
            p_all(:,:) = 0.d0
        end if
        allocate(un(0:ny+1,0:nx+1))
        allocate(vn(0:ny+1,0:nx+1))
        allocate(pn(0:ny+1,0:nx+1))
        allocate(b(ny,nx))
        allocate(u(ny,nx))
        allocate(v(ny,nx))
        allocate(p(ny,nx))
        un(:,:) = 0.d0
        vn(:,:) = 0.d0
        pn(:,:) = 0.d0
        b(:,:)  = 0.d0
        u(:,:)  = 0.d0
        v(:,:)  = 0.d0
        p(:,:)  = 0.d0
    end subroutine
    subroutine end_Variables()
        implicit none
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
end module
    
