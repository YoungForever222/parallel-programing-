Program Jordan
    Use aparallel
    Use aVariables
    implicit none
    call init_parallel()
    call init_Variables()
    if (root)write(*,*) "init_Variables() over"
    call Scatter()
    if (root)write(*,*) "Scatter() over"
    call Jordan_Method_MPI()
    if (root)write(*,*) "Jacobi_Method_MPI() over"
    call Output()
    call end_parallel()
end
subroutine get_max(lmax,maxvalue,maxrow,maxcolumn,maxrank)
    Use aparallel
    implicit none
    real(8) :: lmax(0:mpisize-1,0:3),maxvalue,l(0:mpisize-1)
    integer :: maxrow,maxcolumn,maxrank,pos(1)
    l = lmax(:,0)
    pos = MAXLOC(l)
    maxvalue = lmax(pos(1),0)
    maxrow = lmax(pos(1),1)
    maxcolumn = lmax(pos(1),2)
    maxrank = lmax(pos(1),3)
    if (pos(1) .ne. maxrank) then
        write(*,*) "pos .ne. maxrank, MPI_Allgather Error"
        call end_parallel()
        stop
    end if
end subroutine
subroutine Jordan_Method_MPI()
    Use aparallel
    Use aVariables
    implicit none
    integer :: i,j,k,t,v,w,maxrow,maxcolumn,maxrank,shift(n),tem
    real(8) :: lmax(0:mpisize-1,0:3),maxvalue,temp,f(n+1),temp_row(n)
    do i = 1,n
        shift(i) = i
    end do
    do i = 1,m
        do j = 0,mpisize-1
            if(mpirank .lt. j) then     !(1)
                v = i * mpisize + j     !(1.1)
                lmax(mpirank,0) = A_sub(i+1,v)  !(1.2)
                do k = i+1,m            !(1.3)
                    do t = v,n
                        if(abs(A_sub(k,t)) .gt. lmax(mpirank,0)) then
                            lmax(mpirank,0) = A_sub(k,t)
                            lmax(mpirank,1) = k
                            lmax(mpirank,2) = t
                            lmax(mpirank,3) = mpirank
                        end if
                    end do
                end do
            end if
            if(mpirank .ge. j) then     !(2)
                v = i * mpisize + j     !(2.1)
                lmax(mpirank,0) = A_sub(i,v)!(2.2)
                do k = i,m              !(2.3)
                    do t = v,n
                        if(abs(A_sub(k,t)) .gt. lmax(mpirank,0)) then
                            lmax(mpirank,0) = A_sub(k,t)
                            lmax(mpirank,1) = k
                            lmax(mpirank,2) = t
                            lmax(mpirank,3) = mpirank
                        end if
                    end do
                end do
            end if
            call MPI_Barrier(MPI_COMM_WORLD,mpicode)
            call MPI_Allgather(lmax(mpirank,0:3),4,MPI_DP,lmax,4,MPI_DP,MPI_COMM_WORLD,mpicode) !(3)
            call get_max(lmax,maxvalue,maxrow,maxcolumn,maxrank)                                !(4)
            !   (5) 列交换
            if (maxcolumn .ne. v) then
                do t = 1,m                          !(5.1)
                    temp = A_sub(t,v)
                    A_sub(t,v) = A_sub(t,maxcolumn)
                    A_sub(t,maxcolumn) = temp
                end do
                tem = shift(v)                      !(5.2)
                shift(v) = shift(maxcolumn)
                shift(maxcolumn) = tem
            end if
            !   (6) 行交换
            if (mpirank .eq. j) then
                if(maxcolumn .ne. v) then           !书中的判断条件似乎有些问题
                    if(maxrank .eq. j .and. i .ne. maxrow) then 
                        !innerexchangerow()
                        temp_row(:) = A_sub(v,:)
                        A_sub(v,:)  = A_sub(maxrow,:)
                        A_sub(maxrow,:) = temp_row(:)
                        temp = b_sub(v)
                        b_sub(v) = b_sub(maxrow)
                        b_sub(maxrow) =  temp
                    end if
                    if(maxrank .ne. j) then
                        !outerexchangerow()
                        call MPI_Sendrecv(A_sub(v,:),n,MPI_DP,maxrank,mpirank+mpisize, &
                            temp_row(:),n,MPI_DP,maxrank,mpirank+mpisize,MPI_COMM_WORLD,status,mpicode)
                        call MPI_Sendrecv(b(v),1,MPI_DP,maxrank,mpirank+2*mpisize, &
                            temp,1,MPI_DP,maxrank,mpirank+2*mpisize,MPI_COMM_WORLD,status,mpicode)
                        A_sub(v,:) = temp_row(:)
                        b_sub(v) = temp
                    end if
                end if
            end if
            if(mpirank .eq. maxrank .and. maxrank .ne. j) then
                call MPI_Sendrecv(A_sub(maxrow,:),n,MPI_DP,j,j+mpisize, &
                            temp_row(:),n,MPI_DP,j,j+mpisize,MPI_COMM_WORLD,status,mpicode)
                call MPI_Sendrecv(b(maxrow),1,MPI_DP,maxrank,mpirank+2*mpisize, &
                            temp,1,MPI_DP,maxrank,mpirank+2*mpisize,MPI_COMM_WORLD,status,mpicode)
                A_sub(maxrow,:) = temp_row(:)
                b_sub(maxrow) = temp
            end if
            if (mpirank .eq. j) then
                !   (7)
                do k = v+1, n                       !(7.1)
                    A_sub(i,k) = A_sub(i,k)/A_sub(i,v)
                    f(k) = A_sub(i,k)               !(7.3)
                end do
                b_sub(i) = b_sub(i)/A_sub(i,v)      !(7.2)
                A_sub(i,v) = 1.d0                   
                f(n+1) = b_sub(i)                   !(7.4)
                call MPI_Bcast(f,n+1,MPI_DP,j,MPI_COMM_WORLD,mpicode)   !(7.5)
                do k = 1,m
                    if(k .ne. i) then
                        A_sub(k,v) = 0.d0
                        do w = v+1,n                !(i)
                            A_sub(k,w) = A_sub(k,w) - f(w)*A_sub(k,v)
                        end do
                        b_sub(k) = b_sub(k)-f(n+1)*A_sub(k,v)           !(ii)
                    end if
                end do
            else                                    !非主行所在处理器
                call MPI_Bcast(f,n+1,MPI_DP,j,MPI_COMM_WORLD,mpicode)   !(7.7)
                do k = 1,m
                    A_sub(k,v) = 0.d0
                    do w = v+1,n                !(i)
                        A_sub(k,w) = A_sub(k,w) - f(w)*A_sub(k,v)
                    end do
                    b_sub(k) = b_sub(k)-f(n+1)*A_sub(k,v)           !(ii)
                end do
            end if
        end do
    end do
    call MPI_Gather(b_sub,m,MPI_DP,x_1,m,MPI_DP,0,MPI_COMM_WORLD,mpicode)
    if(root) then
        do k = 1,n
            do i = 1,n
                if(shift(i) .eq. k) x(k) = x_1(i)
            end do
        end do
    end if
end subroutine

