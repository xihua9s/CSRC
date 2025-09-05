!> 20221220
!> xihua
!> initial sod problem

subroutine initial_sod_rectangle_ca(supt,den,pre,sound,enin,eos, &
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
	
   supt = 0.2d0 
   gamar = 7.d0 / 5.d0 

	xl=0.d0 ; xr=1.d0
	yd=0.d0 ; yu=0.1d0 
	dx=(xr-xl)/nx  
	dy=(yu-yd)/ny
	
	do  i=1,nx 
	do  j=1,ny   
		if(i <= nx/2) then  
			den(i,j) = 1.d0
			eos(i,j) = 715
			pre(i,j) = 1.d0
		else 
			den(i,j) = 0.125d0
			eos(i,j) = 715
			pre(i,j) = 0.1d0
		endif 
		sound(i,j)= dsqrt( gamar*pre(i,j)/den(i,j) )
		enin(i,j) = pre(i,j)/(gamar-1.d0)/den(i,j)      
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
      bnd_type(i,j) = 100
      Vertex_Type_norm_pre(i,j,1) = 1.d0 
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

end subroutine initial_sod_rectangle_ca
	



subroutine initial_sod_circle_cy_aw(supt,den,pre,sound,enin,eos, &
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
	
   supt = 0.2d0 
   gamar = 7.d0 / 5.d0 

	xl=0.d0 ; xr=1.d0
	yd=0.d0 ; yu=3.1415926535898d0/2.d0 
	dx=(xr-xl)/nx  
	dy=(yu-yd)/ny
	
	do  i=1,nx 
	do  j=1,ny   
		if(i <= nx/2) then  
			den(i,j) = 1.d0
			eos(i,j) = 715
			pre(i,j) = 1.d0
		else 
			den(i,j) = 0.125d0
			eos(i,j) = 715
			pre(i,j) = 0.1d0
		endif 
		sound(i,j)= dsqrt( gamar*pre(i,j)/den(i,j) )
		enin(i,j) = pre(i,j)/(gamar-1.d0)/den(i,j)
	enddo
	enddo 

			
	do i=0,nx
	do j=0,ny
		vertx(i,j)= (xl + dx * i) * cos(yd + dy * j)
		verty(i,j)= (xl + dx * i) * sin(yd + dy * j)
      if ( j == 0 ) then 
         verty(i,j) = 0.d0 
      elseif ( j == ny ) then 
         vertx(i,j) = 0.d0
      elseif ( i == 0) then 
         vertx(i,j) = 0.d0
         verty(i,j) = 0.d0         
      endif  
		vstarx(i,j) = 0.d0 
		vstary(i,j) = 0.d0 
	enddo
	enddo

	do j=1,ny
	do i=1,nx
		area(i,j) = 0.d0
      volm(i,j) = 0.d0

		xr1=vertx(i-1,j-1)           
		yr1=verty(i-1,j-1)
		xr2=vertx(i,j-1)           
		yr2=verty(i,j-1)	            
		area(i,j)=area(i,j)+ca_area(xr1,yr1,xr2,yr2)
		volm(i,j)=volm(i,j)+cy_volm(xr1,yr1,xr2,yr2)

		xr1=vertx(i,j-1)           
		yr1=verty(i,j-1)
		xr2=vertx(i,j)           
		yr2=verty(i,j)	            
		area(i,j)=area(i,j)+ca_area(xr1,yr1,xr2,yr2)
		volm(i,j)=volm(i,j)+cy_volm(xr1,yr1,xr2,yr2)

		xr1=vertx(i,j)           
		yr1=verty(i,j)
		xr2=vertx(i-1,j)           
		yr2=verty(i-1,j)	            
		area(i,j)=area(i,j)+ca_area(xr1,yr1,xr2,yr2)
		volm(i,j)=volm(i,j)+cy_volm(xr1,yr1,xr2,yr2)

		xr1=vertx(i-1,j)           
		yr1=verty(i-1,j)
		xr2=vertx(i-1,j-1)           
		yr2=verty(i-1,j-1)	            
		area(i,j)=area(i,j)+ca_area(xr1,yr1,xr2,yr2)
		volm(i,j)=volm(i,j)+cy_volm(xr1,yr1,xr2,yr2)

      !> 计算子网格的柱几何体积和柱几何质量
      ! call cal_subcell_volm_cy(vertx(i-1,j-1),verty(i-1,j-1), &
      !       vertx(i,j-1),verty(i,j-1),vertx(i,j),verty(i,j),  & 
      !       vertx(i-1,j),verty(i-1,j),subcell_volm_cy(i,j,:))  
      ! smass(i,j,:) = den(i,j) * subcell_volm_cy(i,j,:)

      !> 计算子网格的平几何体积和等效质量
      call cal_subcell_area_ca(vertx(i-1,j-1),verty(i-1,j-1), &
            vertx(i,j-1),verty(i,j-1),vertx(i,j),verty(i,j),  & 
            vertx(i-1,j),verty(i-1,j),subcell_area_ca(i,j,:))  
      subcell_area_cy_aw(i,j,:) = 2.d0/3.d0*subcell_area_ca(i,j,:) + &
                                  1.d0/12.d0*area(i,j)
      subcell_volm_cy_aw(i,j,1) = subcell_area_cy_aw(i,j,1) * verty(i-1,j-1)
      subcell_volm_cy_aw(i,j,2) = subcell_area_cy_aw(i,j,2) * verty(i  ,j-1)
      subcell_volm_cy_aw(i,j,3) = subcell_area_cy_aw(i,j,3) * verty(i  ,j  )
      subcell_volm_cy_aw(i,j,4) = subcell_area_cy_aw(i,j,4) * verty(i-1,j  )
      smass(i,j,:) = den(i,j) * subcell_volm_cy_aw(i,j,:)

      amass(i,j)= sum(smass(i,j,:))
	enddo
	enddo

   do i = 1 , nx - 1
   do j = 1 , ny - 1 
      pmass(i,j) = smass(i,j,3) + smass(i+1,j,4) &
                 + smass(i+1,j+1,1) + smass(i,j+1,2)
      ! pvolm(i,j) = subcell_volm_cy(i,j,3)     + subcell_volm_cy(i+1,j,4) &
      !            + subcell_volm_cy(i+1,j+1,1) + subcell_volm_cy(i,j+1,2) 
      ! parea(i,j) = subcell_area_ca(i,j,3)     + subcell_area_ca(i+1,j,4) &
      !            + subcell_area_ca(i+1,j+1,1) + subcell_area_ca(i,j+1,2) 
      pvolm(i,j) = subcell_volm_cy_aw(i,j,3) + subcell_volm_cy_aw(i+1,j,4) &
                 + subcell_volm_cy_aw(i+1,j+1,1) + subcell_volm_cy_aw(i,j+1,2) 
      parea(i,j) = subcell_area_cy_aw(i,j,3) + subcell_area_cy_aw(i+1,j,4) &
                 + subcell_area_cy_aw(i+1,j+1,1) + subcell_area_cy_aw(i,j+1,2) 
   enddo 
   enddo 

   do i = 0 , 0
   do j = 1 , ny - 1 
      pmass(i,j) = smass(i+1,j,4) + smass(i+1,j+1,1) 
      ! pvolm(i,j) = subcell_volm_cy(i+1,j,4) + subcell_volm_cy(i+1,j+1,1)  
      ! parea(i,j) = subcell_area_ca(i+1,j,4) + subcell_area_ca(i+1,j+1,1)
      pvolm(i,j) = subcell_volm_cy_aw(i+1,j,4) + subcell_volm_cy_aw(i+1,j+1,1) 
      parea(i,j) = subcell_area_cy_aw(i+1,j,4) + subcell_area_cy_aw(i+1,j+1,1) 
   enddo 
   enddo 

   do i = nx , nx 
   do j = 1 , ny - 1 
      pmass(i,j) = smass(i,j,3) + smass(i,j+1,2)
      ! pvolm(i,j) = subcell_volm_cy(i,j,3)     + subcell_volm_cy(i,j+1,2) 
      ! parea(i,j) = subcell_area_ca(i,j,3)     + subcell_area_ca(i,j+1,2) 
      pvolm(i,j) = subcell_volm_cy_aw(i,j,3) + subcell_volm_cy_aw(i,j+1,2)
      parea(i,j) = subcell_area_cy_aw(i,j,3) + subcell_area_cy_aw(i,j+1,2)
   enddo 
   enddo 

   do i = 1 , nx - 1
   do j = 0 , 0 
      pmass(i,j) = smass(i+1,j+1,1) + smass(i,j+1,2)
      ! pvolm(i,j) = subcell_volm_cy(i+1,j+1,1) + subcell_volm_cy(i,j+1,2) 
      ! parea(i,j) = subcell_area_ca(i+1,j+1,1) + subcell_area_ca(i,j+1,2) 
      pvolm(i,j) = subcell_volm_cy_aw(i+1,j+1,1) + subcell_volm_cy_aw(i,j+1,2)
      parea(i,j) = subcell_area_cy_aw(i+1,j+1,1) + subcell_area_cy_aw(i,j+1,2)
   enddo 
   enddo 

   do i = 1 , nx - 1
   do j = ny , ny 
      pmass(i,j) = smass(i,j,3) + smass(i+1,j,4)  
      ! pvolm(i,j) = subcell_volm_cy(i,j,3)     + subcell_volm_cy(i+1,j,4) 
      ! parea(i,j) = subcell_area_ca(i,j,3)     + subcell_area_ca(i+1,j,4)
      pvolm(i,j) = subcell_volm_cy_aw(i,j,3) + subcell_volm_cy_aw(i+1,j,4)                 
      parea(i,j) = subcell_area_cy_aw(i,j,3) + subcell_area_cy_aw(i+1,j,4)                 
   enddo 
   enddo 

   do i = 0 , 0
   do j = 0 , 0  
      pmass(i,j) = smass(i+1,j+1,1)
      ! pvolm(i,j) = subcell_volm_cy(i+1,j+1,1) 
      ! parea(i,j) = subcell_area_ca(i+1,j+1,1)  
      pvolm(i,j) = subcell_volm_cy_aw(i+1,j+1,1)
      parea(i,j) = subcell_area_cy_aw(i+1,j+1,1)
   enddo 
   enddo 

   do i = 0,0
   do j = ny,ny 
      pmass(i,j) = smass(i+1,j,4) 
      ! pvolm(i,j) = subcell_volm_cy(i+1,j,4) 
      ! parea(i,j) = subcell_area_ca(i+1,j,4) 
      pvolm(i,j) = subcell_volm_cy_aw(i+1,j,4) 
      parea(i,j) = subcell_area_cy_aw(i+1,j,4) 
   enddo 
   enddo 

   do i = nx,nx
   do j = 0,0
      pmass(i,j) = smass(i,j+1,2)
      ! pvolm(i,j) = subcell_volm_cy(i,j+1,2) 
      ! parea(i,j) = subcell_area_ca(i,j+1,2) 
      pvolm(i,j) = subcell_volm_cy_aw(i,j+1,2)
      parea(i,j) = subcell_area_cy_aw(i,j+1,2)
   enddo 
   enddo 

   do i = nx,nx 
   do j = ny,ny  
      pmass(i,j) = smass(i,j,3)
      ! pvolm(i,j) = subcell_volm_cy(i,j,3) 
      ! parea(i,j) = subcell_area_ca(i,j,3)
      pvolm(i,j) = subcell_volm_cy_aw(i,j,3)
      parea(i,j) = subcell_area_cy_aw(i,j,3)
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
      bnd_type(i,j) = 100
      Vertex_Type_norm_pre(i,j,1) = 0.d0 
      Vertex_Type_norm_pre(i,j,2) = 1.d0 
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

end subroutine initial_sod_circle_cy_aw
	
