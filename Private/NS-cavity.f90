program cavity_flow
    USE aparallel
    USE aVariables
    implicit none
    
    call init_parallel()
    call check_parallel()
    call init_Variables()
    call field_output()
    do n = 1,nt
        if(root) write(*,*) n
        call velocity_interface()
        call build_up_b()
        call pressure_poisson()
        call velocity_update()
        !write(*,*) "mpirank=",mpirank,"un=",un
        if (n/ntout*ntout .eq. n) then
            call field_output()
            !write(*,*) "field_output()"
        end if
    end do
    call end_Variables()
    call end_parallel()
end program





    
