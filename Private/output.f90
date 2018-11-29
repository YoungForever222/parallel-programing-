subroutine field_output()
    USE aparallel
    USE aVariables
    implicit none
    integer :: i,j,rank,x,y
    character*10:: fileName
    integer:: nv
    integer,parameter:: namLen=40,idfile=100,numVar=6 
    real(4),parameter:: ZONEMARKER=299.0,EOHMARKER =357.0
    character*40:: ZoneName='ZONE 1',title="Binary File.", &
            varname(numVar)=(/'x','y','p','u','v','vort'/)
    !   Gather
    !write(*,*) "Gather Start"
    if(root) then
        do i = 0,p_sqrt-1
            do j = 0,p_sqrt-1
                rank = j*p_sqrt+i
                if (i.eq.0 .and. j.eq.0) then
                    u_all(j*ny+1:(j+1)*ny,i*nx+1:(i+1)*nx) = u
                    v_all(j*ny+1:(j+1)*ny,i*nx+1:(i+1)*nx) = v
                    p_all(j*ny+1:(j+1)*ny,i*nx+1:(i+1)*nx) = p
                else
                    call MPI_RECV(u_all(j*ny+1:(j+1)*ny,i*nx+1:(i+1)*nx),ny*nx,MPI_DP,rank,rank,MPI_COMM_WORLD,status,mpicode)
                    call MPI_RECV(v_all(j*ny+1:(j+1)*ny,i*nx+1:(i+1)*nx),ny*nx,MPI_DP,rank,rank+mpisize,MPI_COMM_WORLD,status,mpicode)
                    call MPI_RECV(p_all(j*ny+1:(j+1)*ny,i*nx+1:(i+1)*nx),ny*nx,MPI_DP,rank,rank+2*mpisize,MPI_COMM_WORLD,status,mpicode)
                end if
            end do
        end do
    else
        call MPI_SEND(u,ny*nx,MPI_DP,0,mpirank,MPI_COMM_WORLD,mpicode)
        call MPI_SEND(v,ny*nx,MPI_DP,0,mpirank+mpisize,MPI_COMM_WORLD,mpicode)
        call MPI_SEND(p,ny*nx,MPI_DP,0,mpirank+2*mpisize,MPI_COMM_WORLD,mpicode)
    end if
    call MPI_Barrier(MPI_COMM_WORLD,mpicode)
    !write(*,*) "Gather Over"
    !   Output
    if(root) then
        !   File
        write(fileName,'(I6)') n
        fileName = adjustr(fileName)
        do  i=1,10
            if(fileName(i:i)==' ')fileName(i:i)='0'
        enddo
        open(idfile,file='./DatFlow/Flow'//trim(fileName)//'.plt',form='BINARY')
        !   Vorticity
        DO  x=2, xDim-1
	        DO  y=2, yDim-1
	            vor(y,x)=(v_all(y,x+1)-v_all(y,x-1))/(2.0*dx)-(u_all(y+1,x)-u_all(y-1,x))/(2.0*dy)
	        END DO
	    END DO
	    DO	y=2, yDim-1
		    vor(y,1)=2*vor(y,2)-vor(y,3)
		    vor(y,xDim)=2*vor(y,xDim-1)-vor(y,xDim-2)
        END DO
	    DO	x=1,xDim
		    vor(1,x)=2*vor(2,x)-vor(3,x)
		    vor(yDim,x)=2*vor(yDim-1,x)-vor(yDim-2,x)
        END DO
        !   Tecplot File Type
        write(idfile) "#!TDV101"    
        write(idfile) 1
        call dumpstring(title,idfile)
        write(idfile) numVar
        do  nv=1,numVar
            call dumpstring(varname(nv),idfile)
        enddo
        write(idfile) ZONEMARKER	 		 
        call dumpstring(zonename,idfile)
        write(idfile) -1,0,1,0,0,ydim,xDim,1,0
        write(idfile) EOHMARKER
        write(idfile) ZONEMARKER
        do  nv=1,numVar
            write(idfile) 1                                 
        enddo
        write(idfile) 0,-1
        do  x=1, xDim
            do  y=1, yDim
                write(idfile) real(x*dx),real(y*dy),real(p_all(y,x)), &
                            real(u_all(y,x)),real(v_all(y,x)),real(vor(y,x))
            enddo
        enddo
    end if
end subroutine
subroutine dumpstring(instring,idfile)
    implicit none
    character(40) instring
    integer:: nascii,ii,len_temp,idfile
    len_temp=LEN_TRIM(instring)
    do	ii=1,len_temp
	    nascii=ICHAR(instring(ii:ii))
		write(idfile) nascii
    enddo
    write(idfile) 0
    return
end subroutine
