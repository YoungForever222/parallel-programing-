Program main
    USE aparallel
    USE aVariables
    implicit none
    call init_parallel()
    if(root) then
        call next()
        call period_analysis()
    end if