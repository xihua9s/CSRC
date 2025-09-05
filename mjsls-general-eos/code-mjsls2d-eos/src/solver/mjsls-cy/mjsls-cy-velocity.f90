
!> 20230112
!> xihua
!> compute velocity and acceleration based on the force 
   
subroutine cal_node_velocity_cy(nx,ny,mass_p,dt,bnd_type,Vertex_Type_norm_pre, &
            vstarx,vstary,vstarx_pre,vstary_pre,fx_pre,fy_pre,g,vertx,verty,smass)
            
   implicit none 
   integer   :: nx,ny
   real*8    :: vstarx(0:nx,0:ny),vstary(0:nx,0:ny),mass_p(0:nx,0:ny)
   real*8    :: vstarx_pre(0:nx,0:ny),vstary_pre(0:nx,0:ny)
   real*8    :: fx_pre(0:nx,0:ny),fy_pre(0:nx,0:ny),dt,g,smass(1:nx,1:ny,1:4),smass_tri(1:nx,1:ny,1:8)
	integer   :: bnd_type(0:nx,0:ny)
	real*8    :: Vertex_Type_norm_pre(0:nx,0:ny,3)   
	real*8 :: vertx(0:nx,0:ny),verty(0:nx,0:ny)
   
   integer  :: i,j
   real*8   :: massp,ca_area,cy_volm,xr1,yr1,xr2,yr2
   real*8   :: area(1:nx,1:ny),volm(1:nx,1:ny),pvolm(0:nx,0:ny)
   real*8  :: subcell_area_ca(1:nx,1:ny,1:4),parea(0:nx,0:ny)  
   real*8  :: subcell_area_cy_aw(1:nx,1:ny,1:4)    
   real*8  :: subcell_volm_cy_aw(1:nx,1:ny,1:4)    



   vstarx_pre = 0.d0 
   vstary_pre = 0.d0 
   
   !> [110,slide_plane] Slide on a plane
   !> [100,slide_line] Slide on a line
   !> [90,given_pre] Given pressure 1.d-6   
   !> [140,given_velocity] Given velocity
   !> [250,mix_slide_plane_given_pre] Mixed, slide on a plane + Given pressure   
   !> [240,241,242,mix_slide_line_given_pre] Mixed, slide on a line + Given pressure   
   !> [30,fix] Fixed
   !> [50,inner] inner node   
   do j = 1 , ny 
   do i = 0 , nx 
      select case (bnd_type(i,j))
      case (50)
         vstarx_pre(i,j) = vstarx(i,j) - dt * fx_pre(i,j)* verty(i,j)/mass_p(i,j)
         vstary_pre(i,j) = vstary(i,j) - dt * fy_pre(i,j)* verty(i,j)/mass_p(i,j)
      case (90)
         vstarx_pre(i,j) = vstarx(i,j) - dt * fx_pre(i,j)* verty(i,j)/mass_p(i,j)
         vstary_pre(i,j) = vstary(i,j) - dt * fy_pre(i,j)* verty(i,j)/mass_p(i,j)
      case (100)      
         vstarx_pre(i,j) = vstarx(i,j) - dt * Vertex_Type_norm_pre(i,j,1)/mass_p(i,j)* &
         (fx_pre(i,j)*Vertex_Type_norm_pre(i,j,1) + fy_pre(i,j)*Vertex_Type_norm_pre(i,j,2))* verty(i,j)
         vstary_pre(i,j) = vstary(i,j) - dt * Vertex_Type_norm_pre(i,j,2)/mass_p(i,j)* &
         (fx_pre(i,j)*Vertex_Type_norm_pre(i,j,1) + fy_pre(i,j)*Vertex_Type_norm_pre(i,j,2))* verty(i,j)  
      case (30)      
         vstarx_pre(i,j) = vstarx(i,j) - dt * 0.d0 
         vstary_pre(i,j) = vstary(i,j) - dt * 0.d0 
      case (140)      
         vstarx_pre(i,j) = Vertex_Type_norm_pre(i,j,1)
         vstary_pre(i,j) = Vertex_Type_norm_pre(i,j,2) 
      case (240)      
         vstarx_pre(i,j) = vstarx(i,j) - dt * Vertex_Type_norm_pre(i,j,1)/mass_p(i,j)* &
         (fx_pre(i,j)*Vertex_Type_norm_pre(i,j,1) + fy_pre(i,j)*Vertex_Type_norm_pre(i,j,2))* verty(i,j)
         vstary_pre(i,j) = vstary(i,j) - dt * Vertex_Type_norm_pre(i,j,2)/mass_p(i,j)* &
         (fx_pre(i,j)*Vertex_Type_norm_pre(i,j,1) + fy_pre(i,j)*Vertex_Type_norm_pre(i,j,2))* verty(i,j)     
      case (241)      
         vstarx_pre(i,j) = vstarx(i,j) - dt * Vertex_Type_norm_pre(i,j,1)/mass_p(i,j)* &
         (fx_pre(i,j)*Vertex_Type_norm_pre(i,j,1) + fy_pre(i,j)*Vertex_Type_norm_pre(i,j,2))* verty(i,j)
         vstary_pre(i,j) = vstary(i,j) - dt * Vertex_Type_norm_pre(i,j,2)/mass_p(i,j)* &
         (fx_pre(i,j)*Vertex_Type_norm_pre(i,j,1) + fy_pre(i,j)*Vertex_Type_norm_pre(i,j,2))* verty(i,j)     
      case (242)      
         vstarx_pre(i,j) = vstarx(i,j) - dt * Vertex_Type_norm_pre(i,j,1)/mass_p(i,j)* &
         (fx_pre(i,j)*Vertex_Type_norm_pre(i,j,1) + fy_pre(i,j)*Vertex_Type_norm_pre(i,j,2))* verty(i,j)
         vstary_pre(i,j) = vstary(i,j) - dt * Vertex_Type_norm_pre(i,j,2)/mass_p(i,j)* &
         (fx_pre(i,j)*Vertex_Type_norm_pre(i,j,1) + fy_pre(i,j)*Vertex_Type_norm_pre(i,j,2))* verty(i,j)     
      case default
         print*, bnd_type(i,j)
         stop 'bnd_type(i,j) is errorfea'
      end select
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

      ! smass(i,j,:) = den(i,j) * subcell_volm_cy(i,j,:)    !> 子网格真实质量
      ! smass(i,j,:) = den(i,j) * subcell_volm_cy_aw(i,j,:) !> 子网格面质量
      ! amass(i,j)= sum(smass(i,j,:))

	enddo
	enddo

   do i = 1 , nx - 1
   do j = 1 , ny - 1 
      pvolm(i,j) = subcell_volm_cy_aw(i,j,3) + subcell_volm_cy_aw(i+1,j,4) &
                 + subcell_volm_cy_aw(i+1,j+1,1) + subcell_volm_cy_aw(i,j+1,2) 
      parea(i,j) = subcell_area_cy_aw(i,j,3) + subcell_area_cy_aw(i+1,j,4) &
                 + subcell_area_cy_aw(i+1,j+1,1) + subcell_area_cy_aw(i,j+1,2) 
   enddo 
   enddo 

   do i = 0 , 0
   do j = 1 , ny - 1 
      pvolm(i,j) = subcell_volm_cy_aw(i+1,j,4) + subcell_volm_cy_aw(i+1,j+1,1) 
      parea(i,j) = subcell_area_cy_aw(i+1,j,4) + subcell_area_cy_aw(i+1,j+1,1) 
   enddo 
   enddo 

   do i = nx , nx 
   do j = 1 , ny - 1 
      pvolm(i,j) = subcell_volm_cy_aw(i,j,3) + subcell_volm_cy_aw(i,j+1,2)
      parea(i,j) = subcell_area_cy_aw(i,j,3) + subcell_area_cy_aw(i,j+1,2)
   enddo 
   enddo 

   do i = 1 , nx - 1
   do j = 0 , 0 
      pvolm(i,j) = subcell_volm_cy_aw(i+1,j+1,1) + subcell_volm_cy_aw(i,j+1,2)
      parea(i,j) = subcell_area_cy_aw(i+1,j+1,1) + subcell_area_cy_aw(i,j+1,2)
   enddo 
   enddo 

   do i = 1 , nx - 1
   do j = ny , ny               
      pvolm(i,j) = subcell_volm_cy_aw(i,j,3) + subcell_volm_cy_aw(i+1,j,4)                 
      parea(i,j) = subcell_area_cy_aw(i,j,3) + subcell_area_cy_aw(i+1,j,4)                 
   enddo 
   enddo 

   do i = 0 , 0
   do j = 0 , 0  
      pvolm(i,j) = subcell_volm_cy_aw(i+1,j+1,1)
      parea(i,j) = subcell_area_cy_aw(i+1,j+1,1)
   enddo 
   enddo 

   do i = 0,0
   do j = ny,ny 
      pvolm(i,j) = subcell_volm_cy_aw(i+1,j,4) 
      parea(i,j) = subcell_area_cy_aw(i+1,j,4) 
   enddo 
   enddo 

   do i = nx,nx
   do j = 0,0
      pvolm(i,j) = subcell_volm_cy_aw(i,j+1,2)
      parea(i,j) = subcell_area_cy_aw(i,j+1,2)
   enddo 
   enddo 

   do j = 0 , 0 
   do i = 0 , nx 
      if( i == 0 ) then 
         massp = mass_p(i,1)/pvolm(i,1)* &
                  (  subcell_area_cy_aw(i+1,j+1,1) )
      elseif ( i == nx ) then 
         massp = mass_p(i,1)/pvolm(i,1) * &
                  ( subcell_area_cy_aw(i,j+1,2) )
      else 
         massp = mass_p(i,1)/pvolm(i,1) * &
                  ( subcell_area_cy_aw(i,j+1,2) + subcell_area_cy_aw(i+1,j+1,1) )
      endif  

      ! if( i == 0 ) then 
      !    massp = smass(i+1,j+1,4)/subcell_volm_cy_aw(i+1,j+1,4)*subcell_area_cy_aw(i+1,j+1,1)
      ! elseif ( i == nx ) then 
      !    massp = smass(i,j+1,3)/subcell_volm_cy_aw(i,j+1,3) * subcell_area_cy_aw(i,j+1,2)
      ! else 
      !    massp = smass(i,j+1,3)/subcell_volm_cy_aw(i,j+1,3) * subcell_area_cy_aw(i,j+1,2) + &
      !            smass(i+1,j+1,4)/subcell_volm_cy_aw(i+1,j+1,4)*subcell_area_cy_aw(i+1,j+1,1)
      ! endif  

      select case (bnd_type(i,j))
      case (50)
         vstarx_pre(i,j) = vstarx(i,j) - dt * fx_pre(i,j)/massp
         vstary_pre(i,j) = vstary(i,j) - dt * fy_pre(i,j)/massp
      case (90)
         vstarx_pre(i,j) = vstarx(i,j) - dt * fx_pre(i,j)/massp
         vstary_pre(i,j) = vstary(i,j) - dt * fy_pre(i,j)/massp
      case (100)      
         vstarx_pre(i,j) = vstarx(i,j) - dt * Vertex_Type_norm_pre(i,j,1)/massp* &
         (fx_pre(i,j)*Vertex_Type_norm_pre(i,j,1) + fy_pre(i,j)*Vertex_Type_norm_pre(i,j,2))
         vstary_pre(i,j) = vstary(i,j) - dt * Vertex_Type_norm_pre(i,j,2)/massp* &
         (fx_pre(i,j)*Vertex_Type_norm_pre(i,j,1) + fy_pre(i,j)*Vertex_Type_norm_pre(i,j,2))   
      case (30)      
         vstarx_pre(i,j) = vstarx(i,j) - dt * 0.d0 
         vstary_pre(i,j) = vstary(i,j) - dt * 0.d0 
      case (140)      
         vstarx_pre(i,j) = Vertex_Type_norm_pre(i,j,1)
         vstary_pre(i,j) = Vertex_Type_norm_pre(i,j,2) 
      case (240)      
         vstarx_pre(i,j) = vstarx(i,j) - dt * Vertex_Type_norm_pre(i,j,1)/massp* &
         (fx_pre(i,j)*Vertex_Type_norm_pre(i,j,1) + fy_pre(i,j)*Vertex_Type_norm_pre(i,j,2))
         vstary_pre(i,j) = vstary(i,j) - dt * Vertex_Type_norm_pre(i,j,2)/massp* &
         (fx_pre(i,j)*Vertex_Type_norm_pre(i,j,1) + fy_pre(i,j)*Vertex_Type_norm_pre(i,j,2))     
      case (241)      
         vstarx_pre(i,j) = vstarx(i,j) - dt * Vertex_Type_norm_pre(i,j,1)/massp* &
         (fx_pre(i,j)*Vertex_Type_norm_pre(i,j,1) + fy_pre(i,j)*Vertex_Type_norm_pre(i,j,2))
         vstary_pre(i,j) = vstary(i,j) - dt * Vertex_Type_norm_pre(i,j,2)/massp* &
         (fx_pre(i,j)*Vertex_Type_norm_pre(i,j,1) + fy_pre(i,j)*Vertex_Type_norm_pre(i,j,2))     
      case (242)      
         vstarx_pre(i,j) = vstarx(i,j) - dt * Vertex_Type_norm_pre(i,j,1)/massp* &
         (fx_pre(i,j)*Vertex_Type_norm_pre(i,j,1) + fy_pre(i,j)*Vertex_Type_norm_pre(i,j,2))
         vstary_pre(i,j) = vstary(i,j) - dt * Vertex_Type_norm_pre(i,j,2)/massp* &
         (fx_pre(i,j)*Vertex_Type_norm_pre(i,j,1) + fy_pre(i,j)*Vertex_Type_norm_pre(i,j,2))     
      case default
         print*, bnd_type(i,j)
         stop 'bnd_type(i,j) is errorfea'
      end select
   enddo
   enddo 

   do j = 1 , ny-1 
   do i = 0 , nx    
      vstary_pre(i,j) = vstary_pre(i,j) - dt * g
   enddo
   enddo 

end subroutine cal_node_velocity_cy
 

