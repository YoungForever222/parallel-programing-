program cavity_flow
    USE aparallel
    USE aVariables
    implicit none
    call init_parallel()
    call check_parallel()
    
    call end_Variables()
    call end_parallel()
end program

subroutine set_boundary()
    USE aVariables
