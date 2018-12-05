module aVariables
    USE apparallel
    implicit none
    integer,parameter :: N = 1000, M=10
    character(len=N):: T
    character(len=M):: P
    integer         :: newnext(1:M)
    integer         :: next(1:M+1)
    integer         :: period_len,period_num,pattern_suffixlen
contains
    subroutine next()
        implicit none
        integer :: i,j
        next(1) = 0                                 !(1)
        newnext(1) = 0                              !(1)
        j=2                                         !(2)
        do while(j .le. M+1)                        !(3)
            i = next(j-1)                           !(3.1)
            do while(i .ne. 0 .and. P(i) .ne. P(j-1))!(3.2)
                i = next(i)
            end do
            next(j) = i+1                           !(3.3)
            if (j .ne. M+1) then                    !(3.4)
                if(P(j) .ne. P(j+1)) then
                    newnext(j) = i+1
                else
                    newnext(j) = newnext(i+1)
                end if
            end if
            j = j + 1                               !(3.5)      
        end do
    subroutine period_analysis()
        implicit none
        period_len = m + 1 - next(m + 1)
        period_num = int(m/period)
        pattern_suffixlen = mod(m,period_len)
        