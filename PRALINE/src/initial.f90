


subroutine initial_sls_triple(supt,den,pre,sound,enin,eos, &
       vertx,verty,amass,gamar,nx,ny,vstarx,vstary,bnd_type,Vertex_Type_norm_pre)
	implicit none
	real*8 :: den(1:nx,1:ny),pre(1:nx,1:ny),amass(1:nx,1:ny)
	real*8 :: sound(1:nx,1:ny),enin(1:nx,1:ny)
	real*8 :: area(1:nx,1:ny)
	real*8 :: vertx(0:nx,0:ny),verty(0:nx,0:ny)
	real*8 :: vstarx(0:nx,0:ny),vstary(0:nx,0:ny)
	integer  :: eos(1:nx,1:ny),bnd_type(0:nx,0:ny)
	real*8  :: Vertex_Type_norm_pre(0:nx,0:ny,3)
	integer :: nx,ny
	real*8  :: gamar,supt
	
	real*8 :: PI,xl,xr,yd,yu,dx,dy,r,theta,xr1,yr1,xr2,yr2,x,y
	integer :: i,j
	
	supt  = 5.d0 
   gamar = 7.d0/5.d0 
   Pi    = 3.14159265359d0 
   
	xl=0.d0 ; xr=7.d0
	yd=0.d0 ; yu=3.d0
	dx=(xr-xl)/nx  
	dy=(yu-yd)/ny
	
	do i=0,nx
	do j=0,ny
		vertx(i,j)= xl + dx * i
		verty(i,j)= yd + dy * j
	enddo
   enddo   
   vstarx = 0.d0 
   vstary = 0.d0
   
	do  i=1,nx 
	do  j=1,ny    
      x = ( vertx(i-1,j-1) + vertx(i-1,j) + vertx(i,j-1) + vertx(i,j) ) / 4.d0
      y = ( verty(i-1,j-1) + verty(i-1,j) + verty(i,j-1) + verty(i,j) ) / 4.d0
      if ( x < 1.d0 ) then 
         den(i,j) = 1.d0
         pre(i,j) = 1.d0
         eos(i,j) = 312 !>1.5
      else 
         if ( y < 1.5 ) then 
            den(i,j) = 1.d0
            pre(i,j) = 0.1d0 
            eos(i,j) = 715 !>1.4
         else 
            den(i,j) = 0.125d0
            pre(i,j) = 0.1d0 
            eos(i,j) = 815 !>1.6
         endif 
      endif 
      call cal_eos(eos(i,j),gamar)
      sound(i,j)= dsqrt( gamar*pre(i,j)/den(i,j) )
      enin(i,j) = pre(i,j)/(gamar-1.d0)/den(i,j)      
	enddo
	enddo 
			
         
	do i=1,nx
	do j=1,ny
		area(i,j)=0.0

		xr1=vertx(i-1,j-1)           
		yr1=verty(i-1,j-1)
		xr2=vertx(i,j-1)           
		yr2=verty(i,j-1)	            
		area(i,j)=area(i,j)+0.5*(xr1*yr2-xr2*yr1)

		xr1=vertx(i,j-1)           
		yr1=verty(i,j-1)
		xr2=vertx(i,j)           
		yr2=verty(i,j)	            
		area(i,j)=area(i,j)+0.5*(xr1*yr2-xr2*yr1)

		xr1=vertx(i,j)           
		yr1=verty(i,j)
		xr2=vertx(i-1,j)           
		yr2=verty(i-1,j)	            
		area(i,j)=area(i,j)+0.5*(xr1*yr2-xr2*yr1)

		xr1=vertx(i-1,j)           
		yr1=verty(i-1,j)
		xr2=vertx(i-1,j-1)           
		yr2=verty(i-1,j-1)	            
		area(i,j)=area(i,j)+0.5*(xr1*yr2-xr2*yr1)

		amass(i,j)=den(i,j)*area(i,j)      
	enddo
	enddo

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

end subroutine initial_sls_triple
	



