!> 20221220
!> xihua
!> initial rti problem
subroutine initial_rti_ca(supt,den,pre,sound,enin,eos, &
   vertx,verty,amass,gamar,nx,ny,vstarx,vstary,bnd_type,&
   Vertex_Type_norm_pre,pmass,smass,g)
	implicit none
	real*8 :: den(1:nx,1:ny),pre(1:nx,1:ny),amass(1:nx,1:ny)
	real*8 :: sound(1:nx,1:ny),enin(1:nx,1:ny)
	real*8 :: area(1:nx,1:ny),volm(1:nx,1:ny)
	real*8 :: vertx(0:nx,0:ny),verty(0:nx,0:ny),pmass(0:nx,0:ny)
	real*8 :: vstarx(0:nx,0:ny),vstary(0:nx,0:ny)
	integer  :: eos(1:nx,1:ny),bnd_type(0:nx,0:ny)
	real*8  :: Vertex_Type_norm_pre(0:nx,0:ny,3)
	integer :: nx,ny
	real*8  :: gamar,supt,ca_area,cy_volm,parea(0:nx,0:ny),pvolm(0:nx,0:ny)
   real*8  :: subcell_volm_cy(1:nx,1:ny,1:4),smass(1:nx,1:ny,1:4),smass_tri(1:nx,1:ny,1:8)
   real*8  :: subcell_area_ca(1:nx,1:ny,1:4) ,g
   real*8  :: subcell_area_cy_aw(1:nx,1:ny,1:4)    
   real*8  :: subcell_volm_cy_aw(1:nx,1:ny,1:4) 
	
	real*8 :: PI,xl,xr,yd,yu,dx,dy,r,theta,xr1,yr1,xr2,yr2,x,y
	integer :: i,j
	
	supt  = 10.d0 
   gamar = 7.d0/5.d0 
   Pi    = 3.14159265359d0 
   g     = 0.1d0  !> gravity 
   
	xl=0.d0 ; xr=1.d0/3.d0
	yd=0.d0 ; yu=1.d0
	! xl=0.d0 ; xr=1.d0
	! yd=0.d0 ; yu=6.d0
	dx=(xr-xl)/nx  
	dy=(yu-yd)/ny

	do i=0,nx
	do j=ny/2 , ny/2
		vertx(i,j)= xl + dx * i
		verty(i,j)= 0.5d0 + 1.d0/100.d0*cos(6.d0*PI*vertx(i,j))
		! verty(i,j)= 3.d0 + 0.06d0 * cos(2.d0*PI*vertx(i,j))
	enddo
   enddo   

	do i=0,nx
	do j=0,ny/2-1
		vertx(i,j)= xl + dx * i
		verty(i,j)= yd + (verty(i,ny/2)-yd)/(ny/2) * j
	enddo
   enddo   
   
	do i=0,nx
	do j=ny/2+1 , ny 
		vertx(i,j)= xl + dx * i
		verty(i,j)= verty(i,ny/2) + (yu - verty(i,ny/2))/(ny/2) * (j-ny/2)
	enddo
   enddo   
   vstarx = 0.d0 
   vstary = 0.d0 

	do  i=1,nx 
	do  j=1,ny    
      x = ( vertx(i-1,j-1) + vertx(i-1,j) + vertx(i,j-1) + vertx(i,j) ) / 4.d0
      y = ( verty(i-1,j-1) + verty(i-1,j) + verty(i,j-1) + verty(i,j) ) / 4.d0
      ! if ( j <= ny/2 ) then 
      if( y <= verty(i,ny/2) ) then 
         den(i,j) = 1.d0
         pre(i,j) = 1.d0 + 2.d0*g*(1.d0-verty(i,ny/2)) + 1.d0*g*(verty(i,ny/2)-y)
         eos(i,j) = 715
      else
         den(i,j) = 2.d0
         pre(i,j) = 1.d0 + 2.d0*g*(1.d0-y)
         eos(i,j) = 715
      endif 
      call cal_eos(eos(i,j),gamar)
      sound(i,j)= dsqrt( gamar*pre(i,j)/den(i,j) )
      enin(i,j) = pre(i,j)/(gamar-1.d0)/den(i,j)      
	enddo
	enddo 

   call cal_geo_2d(nx,ny,vertx,verty,den,area,volm,&
         subcell_area_ca,smass,amass,parea,pvolm,pmass,smass_tri)

   !> Boundary condition:
   !> [110,slide_plane] Slide on a plane
   !> [100,slide_line] Slide on a line
   !> [90,given_pre] Given pressure 1.d-6   
   !> [140,given_velocity] Given velocity
   !> [250,mix_slide_plane_given_pre] Mixed, slide on a plane + Given pressure   
   !> [240,241,242,mix_slide_line_given_pre] Mixed, slide on a line + Given pressure   
   !> [30,fix] Fixed
   !> [50,inner] inner node
   
   bnd_type = 50
   Vertex_Type_norm_pre = 0.d0 
   
   j = 0
   do i = 1 , nx-1
      bnd_type(i,j) = 100
      Vertex_Type_norm_pre(i,j,1) = 1.d0 
      Vertex_Type_norm_pre(i,j,2) = 0.d0 
      Vertex_Type_norm_pre(i,j,3) = 0.d0 
   enddo 
   j = ny
   do i = 1 , nx-1
      bnd_type(i,j) = 100
      Vertex_Type_norm_pre(i,j,1) = 1.d0 
      Vertex_Type_norm_pre(i,j,2) = 0.d0 
      Vertex_Type_norm_pre(i,j,3) = 0.d0 
   enddo 
   i = 0
   do j = 1 , ny-1 
      bnd_type(i,j) = 100
      Vertex_Type_norm_pre(i,j,1) = 0.d0 
      Vertex_Type_norm_pre(i,j,2) = 1.d0 
      Vertex_Type_norm_pre(i,j,3) = 0.d0 
   enddo
   i = nx 
   do j = 1 , ny-1 
      bnd_type(i,j) = 100
      Vertex_Type_norm_pre(i,j,1) = 0.d0 
      Vertex_Type_norm_pre(i,j,2) = 1.d0 
      Vertex_Type_norm_pre(i,j,3) = 0.d0 
   enddo

   i = 0
   j = 0 
   bnd_type(i,j) = 30
   Vertex_Type_norm_pre(i,j,1) = 0.d0 
   Vertex_Type_norm_pre(i,j,2) = 1.d0 
   Vertex_Type_norm_pre(i,j,3) = 0.d0    

   i = nx
   j = 0 
   bnd_type(i,j) = 30
   Vertex_Type_norm_pre(i,j,1) = 0.d0 
   Vertex_Type_norm_pre(i,j,2) = 1.d0 
   Vertex_Type_norm_pre(i,j,3) = 0.d0  
   
   i = 0
   j = ny
   bnd_type(i,j) = 30
   Vertex_Type_norm_pre(i,j,1) = 0.d0 
   Vertex_Type_norm_pre(i,j,2) = 1.d0 
   Vertex_Type_norm_pre(i,j,3) = 0.d0  

   i = nx
   j = ny
   bnd_type(i,j) = 30
   Vertex_Type_norm_pre(i,j,1) = 0.d0 
   Vertex_Type_norm_pre(i,j,2) = 1.d0 
   Vertex_Type_norm_pre(i,j,3) = 0.d0  

	!> Growth
	 !do it = 1 , 1000
		 !t = 6.d0/1000 * it  
		 !write(12,*) t , 0.02d0*exp( sqrt(1.d0/3.d0*0.1d0*1.5d0*3.1416926d0/(2.d0/3.d0))*t )
	 !enddo 
	 !stop
end subroutine initial_rti_ca
	

