

subroutine  output_sls(nx,ny,den,pre,enin, &
            vertx,verty,gamar,vstarx,vstary,it)
	implicit none
	integer,intent(in) :: nx,ny,it
	real*8 :: den(1:nx,1:ny),pre(1:nx,1:ny),enin(1:nx,1:ny)
	real*8 :: vertx(0:nx,0:ny),verty(0:nx,0:ny)
	real*8 :: vstarx(0:nx,0:ny),vstary(0:nx,0:ny)
	real*8  :: gamar
	character*8  :: fname 

	integer  :: i,j 

	write(fname,'(I8.8)') it
	open(18,file='phys_it_'//fname//'.plt',status='unknown') 
	write(18,*) 'title="contour"'
	write(18,*) 'variables="x","y","vx","vy","den","pre","enin"'
	write(18,*)  'zone  i=', nx+1, 'j=' ,ny+1, ',F=BLOCK,VARLOCATION=([5-7]=CELLCENTERED)'
	do j=0,ny
	do i=0,nx
		write(18,'(f20.8)') vertx(i,j)
	enddo 
	enddo
	do j=0,ny
	do i=0,nx
		write(18,'(f20.8)') verty(i,j)
	enddo 
	enddo
	do j=0,ny
	do i=0,nx
		write(18,'(f20.8)') vstarx(i,j)
	enddo 
	enddo
	do j=0,ny
	do i=0,nx
		write(18,'(f20.8)') vstary(i,j)
	enddo 
	enddo
	do j=1,ny
	do i=1,nx
		write(18,'(f20.8)') den(i,j)
	enddo 
	enddo
	do j=1,ny
	do i=1,nx
		write(18,'(f20.8)') pre(i,j)
	enddo 
	enddo
	do j=1,ny
	do i=1,nx
		write(18,'(f20.8)') enin(i,j)
	enddo 
	enddo
	close(18)

end subroutine  output_sls


subroutine  output_sls_final(nx,ny,den,pre,enin, &
            vertx,verty,gamar,vstarx,vstary,it)
	implicit none
	integer,intent(in) :: nx,ny,it
	real*8 :: den(1:nx,1:ny),pre(1:nx,1:ny),enin(1:nx,1:ny)
	real*8 :: vertx(0:nx,0:ny),verty(0:nx,0:ny)
	real*8 :: vstarx(0:nx,0:ny),vstary(0:nx,0:ny)
	real*8  :: gamar,rx,ry
	character*8  :: fname 

	integer  :: i,j 

	write(fname,'(I8.8)') it
	open(18,file='final.plt',status='unknown') 
	write(18,*) 'title="contour"'
	write(18,*) 'variables="x","y","vx","vy","den","pre","enin"'
	write(18,*)  'zone  i=', nx+1, 'j=' ,ny+1, ',F=BLOCK,VARLOCATION=([5-7]=CELLCENTERED)'
	do j=0,ny
	do i=0,nx
		write(18,'(f20.8)') vertx(i,j)
	enddo 
	enddo
	do j=0,ny
	do i=0,nx
		write(18,'(f20.8)') verty(i,j)
	enddo 
	enddo
	do j=0,ny
	do i=0,nx
		write(18,'(f20.8)') vstarx(i,j)
	enddo 
	enddo
	do j=0,ny
	do i=0,nx
		write(18,'(f20.8)') vstary(i,j)
	enddo 
	enddo
	do j=1,ny
	do i=1,nx
		write(18,'(f20.8)') den(i,j)
	enddo 
	enddo
	do j=1,ny
	do i=1,nx
		write(18,'(f20.8)') pre(i,j)
	enddo 
	enddo
	do j=1,ny
	do i=1,nx
		write(18,'(f20.8)') enin(i,j)
	enddo 
	enddo
	close(18)

   
	open(18,file='OutputData/final_R.plt',status='unknown') 
	write(18,*) 'title="contour"'
	write(18,*) 'variables="r","den" '
	write(18,*)  'zone  i=', nx*ny, ',F=POINT'
	do j=1,ny
	do i=1,nx
      rx = ( vertx(i-1,j-1) + vertx(i-0,j-1) + vertx(i-1,j-0) + vertx(i-0,j-0) ) / 4.d0 
      ry = ( verty(i-1,j-1) + verty(i-0,j-1) + verty(i-1,j-0) + verty(i-0,j-0) ) / 4.d0 
		write(18,'(f20.8)') sqrt(rx*rx+ry*ry), den(i,j)
	enddo 
	enddo
	
	close(18)   
   
   
end subroutine  output_sls_final

