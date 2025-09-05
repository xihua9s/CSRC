
  
!> 20230112
!> xihua
!> compute velocity and acceleration based on the force 
   

subroutine cal_node_velocity_ca(nx,ny,mass_p,dt,bnd_type,Vertex_Type_norm_pre, &
            vstarx,vstary,vstarx_pre,vstary_pre,fx_pre,fy_pre,g,vertx,verty,smass)
            
   implicit none 
   integer   :: nx,ny
   real*8    :: vstarx(0:nx,0:ny),vstary(0:nx,0:ny),mass_p(0:nx,0:ny)
   real*8    :: vstarx_pre(0:nx,0:ny),vstary_pre(0:nx,0:ny)
   real*8    :: fx_pre(0:nx,0:ny),fy_pre(0:nx,0:ny),dt,g,smass(1:nx,1:ny,1:4),smass_tri(1:nx,1:ny,1:8)
	integer   :: bnd_type(0:nx,0:ny)
	real*8    :: Vertex_Type_norm_pre(0:nx,0:ny,3)   
	real*8    :: vertx(0:nx,0:ny),verty(0:nx,0:ny)
   
   integer  :: i,j
   real*8   :: massp,ca_area,cy_volm,xr1,yr1,xr2,yr2
   real*8   :: area(1:nx,1:ny),volm(1:nx,1:ny),pvolm(0:nx,0:ny)
   real*8   :: subcell_area_ca(1:nx,1:ny,1:4),parea(0:nx,0:ny)  
   real*8   :: subcell_area_cy_aw(1:nx,1:ny,1:4)    
   real*8   :: subcell_volm_cy_aw(1:nx,1:ny,1:4)    



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
   do j = 0  , ny 
   do i = 0  , nx 
      select case (bnd_type(i,j))
      case (50)
         vstarx_pre(i,j) = vstarx(i,j) - dt * fx_pre(i,j)/mass_p(i,j)
         vstary_pre(i,j) = vstary(i,j) - dt * fy_pre(i,j)/mass_p(i,j)    
      case (90)
         vstarx_pre(i,j) = vstarx(i,j) - dt * fx_pre(i,j)/mass_p(i,j)
         vstary_pre(i,j) = vstary(i,j) - dt * fy_pre(i,j)/mass_p(i,j)
      case (100)      
         vstarx_pre(i,j) = vstarx(i,j) - dt * Vertex_Type_norm_pre(i,j,1)/mass_p(i,j)* &
         (fx_pre(i,j)*Vertex_Type_norm_pre(i,j,1) + fy_pre(i,j)*Vertex_Type_norm_pre(i,j,2))
         vstary_pre(i,j) = vstary(i,j) - dt * Vertex_Type_norm_pre(i,j,2)/mass_p(i,j)* &
         (fx_pre(i,j)*Vertex_Type_norm_pre(i,j,1) + fy_pre(i,j)*Vertex_Type_norm_pre(i,j,2))
      case (30)      
         vstarx_pre(i,j) = vstarx(i,j) - dt * 0.d0 
         vstary_pre(i,j) = vstary(i,j) - dt * 0.d0 
      case (140)      
         vstarx_pre(i,j) = Vertex_Type_norm_pre(i,j,1)
         vstary_pre(i,j) = Vertex_Type_norm_pre(i,j,2) 
      case (240)      
         vstarx_pre(i,j) = vstarx(i,j) - dt * Vertex_Type_norm_pre(i,j,1)/mass_p(i,j)* &
         (fx_pre(i,j)*Vertex_Type_norm_pre(i,j,1) + fy_pre(i,j)*Vertex_Type_norm_pre(i,j,2))
         vstary_pre(i,j) = vstary(i,j) - dt * Vertex_Type_norm_pre(i,j,2)/mass_p(i,j)* &
         (fx_pre(i,j)*Vertex_Type_norm_pre(i,j,1) + fy_pre(i,j)*Vertex_Type_norm_pre(i,j,2))    
      case (241)      
         vstarx_pre(i,j) = vstarx(i,j) - dt * Vertex_Type_norm_pre(i,j,1)/mass_p(i,j)* &
         (fx_pre(i,j)*Vertex_Type_norm_pre(i,j,1) + fy_pre(i,j)*Vertex_Type_norm_pre(i,j,2))
         vstary_pre(i,j) = vstary(i,j) - dt * Vertex_Type_norm_pre(i,j,2)/mass_p(i,j)* &
         (fx_pre(i,j)*Vertex_Type_norm_pre(i,j,1) + fy_pre(i,j)*Vertex_Type_norm_pre(i,j,2))     
      case (242)      
         vstarx_pre(i,j) = vstarx(i,j) - dt * Vertex_Type_norm_pre(i,j,1)/mass_p(i,j)* &
         (fx_pre(i,j)*Vertex_Type_norm_pre(i,j,1) + fy_pre(i,j)*Vertex_Type_norm_pre(i,j,2))
         vstary_pre(i,j) = vstary(i,j) - dt * Vertex_Type_norm_pre(i,j,2)/mass_p(i,j)* &
         (fx_pre(i,j)*Vertex_Type_norm_pre(i,j,1) + fy_pre(i,j)*Vertex_Type_norm_pre(i,j,2))     
      case default
         print*, bnd_type(i,j)
         stop 'bnd_type(i,j) is errorfea'
      end select
      ! print*, i,j,bnd_type(i,j),mass_p(i,j),vstarx_pre(i,j),vstary_pre(i,j) ; pause 121
   enddo
   enddo 

end subroutine cal_node_velocity_ca
 
