!> 20221228
!> xihua
!> initial sedov problem

subroutine initial_sedov_rectangle_ca(supt,den,pre,sound,enin,eos, &
   vertx,verty,amass,gamar,nx,ny,vstarx,vstary,bnd_type,&
   Vertex_Type_norm_pre,pmass,smass,smass_tri)
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
   real*8  :: subcell_area_ca(1:nx,1:ny,1:4)  
   real*8  :: subcell_area_cy_aw(1:nx,1:ny,1:4)    
   real*8  :: subcell_volm_cy_aw(1:nx,1:ny,1:4)   
	
	real*8 :: PI,xl,xr,yd,yu,dx,dy,r,theta,xr1,yr1,xr2,yr2
	integer :: i,j
	
   supt = 1.d0 
   gamar = 7.d0 / 5.d0 

	xl=0.d0 ; xr=1.2d0
	yd=0.d0 ; yu=1.2d0
	dx=(xr-xl)/nx  
	dy=(yu-yd)/ny
	
	do  i=1,nx 
	do  j=1,ny    
		den(i,j) = 1.d0
		pre(i,j) = 1.d-6
		eos(i,j) = 715
	enddo
	enddo 
			
	do i=0,nx
	do j=0,ny
		vertx(i,j)= xl + dx * i
		verty(i,j)= yd + dy * j
		vstarx(i,j) = 0.d0 
		vstary(i,j) = 0.d0 
	enddo
	enddo

   call cal_geo_2d(nx,ny,vertx,verty,den,area,volm,&
         subcell_area_ca,smass,amass,parea,pvolm,pmass,smass_tri)

	enin(1,1) = 0.244816d0 / area(1,1) 
	pre(1,1)  = (gamar-1.d0) * den(1,1) * enin(1,1)   
	do  i=1,nx               
	do  j=1,ny	    
		sound(i,j)= dsqrt( gamar*pre(i,j)/den(i,j) )
		enin(i,j) = pre(i,j)/(gamar-1.d0)/den(i,j)
	enddo
	enddo

	!> Boundary condition:
   !> [110,slide_plane] Slide on a plane
   !> [100,slide_line] Slide on a line
   !> [90,given_pre] Given pressure 1.d-6   
   !> [140,given_velocity] Given velocity
   !> [250,mix_slide_plane_given_pre] Mixed, slide on a plane + Given pressure   
   !> [240,mix_slide_line_given_pre] Mixed, slide on a line + Given pressure   
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
      bnd_type(i,j) = 30
      Vertex_Type_norm_pre(i,j,1) = 0.d0 
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
      bnd_type(i,j) = 30
      Vertex_Type_norm_pre(i,j,1) = 0.d0 
      Vertex_Type_norm_pre(i,j,2) = 0.d0 
      Vertex_Type_norm_pre(i,j,3) = 0.d0 
   enddo

   i = 0
   j = 0 
   bnd_type(i,j) = 30
   Vertex_Type_norm_pre(i,j,1) = 0.d0 
   Vertex_Type_norm_pre(i,j,2) = 0.d0 
   Vertex_Type_norm_pre(i,j,3) = 0.d0    

   i = nx
   j = 0 
   bnd_type(i,j) = 30
   Vertex_Type_norm_pre(i,j,1) = 0.d0 
   Vertex_Type_norm_pre(i,j,2) = 0.d0 
   Vertex_Type_norm_pre(i,j,3) = 0.d0  
   
   i = 0
   j = ny
   bnd_type(i,j) = 30
   Vertex_Type_norm_pre(i,j,1) = 0.d0 
   Vertex_Type_norm_pre(i,j,2) = 0.d0 
   Vertex_Type_norm_pre(i,j,3) = 0.d0  

   i = nx
   j = ny
   bnd_type(i,j) = 30
   Vertex_Type_norm_pre(i,j,1) = 0.d0 
   Vertex_Type_norm_pre(i,j,2) = 0.d0 
   Vertex_Type_norm_pre(i,j,3) = 0.d0     

end subroutine initial_sedov_rectangle_ca
	

!> 20221228
!> xihua
!> initial sedov problem

subroutine initial_sedov4_rectangle_ca(supt,den,pre,sound,enin,eos, &
   vertx,verty,amass,gamar,nx,ny,vstarx,vstary,bnd_type,&
   Vertex_Type_norm_pre,pmass,smass,smass_tri)
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
   real*8  :: subcell_area_ca(1:nx,1:ny,1:4)  
   real*8  :: subcell_area_cy_aw(1:nx,1:ny,1:4)    
   real*8  :: subcell_volm_cy_aw(1:nx,1:ny,1:4) 
	
	real*8 :: PI,xl,xr,yd,yu,dx,dy,r,theta,xr1,yr1,xr2,yr2
	integer :: i,j
	
   supt = 1.d0 
   gamar = 7.d0 / 5.d0 

	xl=-1.2d0 ; xr=1.2d0
	yd=-1.2d0 ; yu=1.2d0
	dx=(xr-xl)/nx  
	dy=(yu-yd)/ny
	
	do  i=1,nx 
	do  j=1,ny    
		den(i,j) = 1.d0
		pre(i,j) = 1.d-6
		eos(i,j) = 715
	enddo
	enddo 
			

	do i=0,nx
	do j=0,ny
      if( i <= nx/5*4 ) then
         dx = (xr-xl)/(nx/5*4)/2.d0 
		   vertx(i,j)= xl + dx * i
      else
         dx = (xr-xl)/(nx/5)/2.d0 
		   vertx(i,j)= 0.d0 +  dx * (i-nx/5*4)
      endif 
      if( j <= ny/5 ) then
         dy = (yu-yd)/(ny/5)/2.d0 
		   verty(i,j)= yd + dy * j
      else
         dy = (yu-yd)/(ny/5*4)/2.d0 
		   verty(i,j)= 0.d0 +  dy * (j-ny/5)
      endif 
		vstarx(i,j) = 0.d0 
		vstary(i,j) = 0.d0 
      ! print*, i,j,verty(i,j)
	enddo
   ! pause 1008
	enddo

	call cal_geo_2d(nx,ny,vertx,verty,den,area,volm,&
         subcell_area_ca,smass,amass,parea,pvolm,pmass,smass_tri)

   j = ny/5
	enin(i,j) = 0.244816d0 / area(i,j) 
	pre(i,j)  = (gamar-1.d0) * den(i,j) * enin(i,j)  

   i = nx/5*4+1 ; j = ny/5
	enin(i,j) = 0.244816d0 / area(i,j) 
	pre(i,j)  = (gamar-1.d0) * den(i,j) * enin(i,j)  

   i = nx/5*4 ; j = ny/5+1
	enin(i,j) = 0.244816d0 / area(i,j) 
	pre(i,j)  = (gamar-1.d0) * den(i,j) * enin(i,j)  

   i = nx/5*4+1 ; j = ny/5+1
	enin(i,j) = 0.244816d0 / area(i,j) 
	pre(i,j)  = (gamar-1.d0) * den(i,j) * enin(i,j)  


	do  i=1,nx               
	do  j=1,ny	    
		sound(i,j)= dsqrt( gamar*pre(i,j)/den(i,j) )
		enin(i,j) = pre(i,j)/(gamar-1.d0)/den(i,j)
	enddo
	enddo

	!> Boundary condition:
   !> [110,slide_plane] Slide on a plane
   !> [100,slide_line] Slide on a line
   !> [90,given_pre] Given pressure 1.d-6   
   !> [140,given_velocity] Given velocity
   !> [250,mix_slide_plane_given_pre] Mixed, slide on a plane + Given pressure   
   !> [240,mix_slide_line_given_pre] Mixed, slide on a line + Given pressure   
   !> [30,fix] Fixed
   !> [50,inner] inner node
   
   bnd_type = 50
   Vertex_Type_norm_pre = 0.d0 
   
   j = 0
   do i = 1 , nx-1
      bnd_type(i,j) = 30
      Vertex_Type_norm_pre(i,j,1) = 0.d0 
      Vertex_Type_norm_pre(i,j,2) = 0.d0 
      Vertex_Type_norm_pre(i,j,3) = 0.d0 
   enddo 
   j = ny
   do i = 1 , nx-1
      bnd_type(i,j) = 30
      Vertex_Type_norm_pre(i,j,1) = 0.d0 
      Vertex_Type_norm_pre(i,j,2) = 0.d0 
      Vertex_Type_norm_pre(i,j,3) = 0.d0 
   enddo 
   i = 0
   do j = 1 , ny-1 
      bnd_type(i,j) = 30
      Vertex_Type_norm_pre(i,j,1) = 0.d0 
      Vertex_Type_norm_pre(i,j,2) = 0.d0 
      Vertex_Type_norm_pre(i,j,3) = 0.d0 
   enddo
   i = nx 
   do j = 1 , ny-1 
      bnd_type(i,j) = 30
      Vertex_Type_norm_pre(i,j,1) = 0.d0 
      Vertex_Type_norm_pre(i,j,2) = 0.d0 
      Vertex_Type_norm_pre(i,j,3) = 0.d0 
   enddo

   i = 0
   j = 0 
   bnd_type(i,j) = 30
   Vertex_Type_norm_pre(i,j,1) = 0.d0 
   Vertex_Type_norm_pre(i,j,2) = 0.d0 
   Vertex_Type_norm_pre(i,j,3) = 0.d0    

   i = nx
   j = 0 
   bnd_type(i,j) = 30
   Vertex_Type_norm_pre(i,j,1) = 0.d0 
   Vertex_Type_norm_pre(i,j,2) = 0.d0 
   Vertex_Type_norm_pre(i,j,3) = 0.d0  
   
   i = 0
   j = ny
   bnd_type(i,j) = 30
   Vertex_Type_norm_pre(i,j,1) = 0.d0 
   Vertex_Type_norm_pre(i,j,2) = 0.d0 
   Vertex_Type_norm_pre(i,j,3) = 0.d0  

   i = nx
   j = ny
   bnd_type(i,j) = 30
   Vertex_Type_norm_pre(i,j,1) = 0.d0 
   Vertex_Type_norm_pre(i,j,2) = 0.d0 
   Vertex_Type_norm_pre(i,j,3) = 0.d0     

end subroutine initial_sedov4_rectangle_ca
	

