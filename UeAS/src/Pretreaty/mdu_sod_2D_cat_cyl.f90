
Module mdu_sod_2D_cat_cyl
   use COMMON_TYPE_EOS
   use COMMON_XXX_Constant
   use COMMON_TYPE_Condition
   use COMMON_TYPE_Mesh
   use COMMON_TYPE_Phys

contains


subroutine Initial_Mesh_Phys_Sod(Ini_CDT,MeshCY2D,PhysCY2D)
   implicit none
   
   Type(COMMON_TYPE_CDT),intent(inout) :: Ini_CDT
   Type(Mesh_Cyl_2D),intent(inout)  :: MeshCY2D
   Type(Phys_Cyl_2D),intent(inout)  :: PhysCY2D
   
   Select Case(Ini_CDT%Gou_Xing)
      
   Case('Sector')
      
      Select Case(Ini_CDT%Angle)

      case('NoEqual')
   
         Call Initial_Mesh_Phys_Sod_Sector_NoEqual_angle_0(Ini_CDT%nx,Ini_CDT%ny,    &
            PhysCY2D%den,PhysCY2D%vel_x,PhysCY2D%vel_y,PhysCY2D%pre,PhysCY2D%sound,  &
            PhysCY2D%tau,PhysCY2D%enin,PhysCY2D%ener,PhysCY2D%EOS,MeshCY2D%vertex_x, &
            MeshCY2D%vertex_y,MeshCY2D%node_type_char,MeshCY2D%node_type_direction, &
            MeshCY2D%node_type_velocity,MeshCY2D%node_type_pressure,Ini_CDT%EOS)
         
      case default
         
         print*, Ini_CDT%Angle
         stop 'Ini_CDT%Angle is Error, pls Check'
         
      End Select
      
   case default

      print*, Ini_CDT%Gou_Xing
      stop 'Ini_CDT%Gou_Xing is Error, pls Check'

   End Select
      


end subroutine Initial_Mesh_Phys_Sod


subroutine Initial_Mesh_Phys_Sod_Sector_NoEqual_angle_0(nx,ny, &
            den,velx,vely,pre,sound,tau,enin,Ener,EOS,&
            vertx,verty,node_type_char,node_type_direction,&
            node_type_velocity,node_type_pressure,eos_num)
   implicit none

   integer,intent(in) :: nx,ny
   integer,intent(in) :: eos_num
   real*8,intent(out) :: den(1:nx,1:ny),velx(1:nx,1:ny),vely(1:nx,1:ny)
   real*8,intent(out) :: pre(1:nx,1:ny),sound(1:nx,1:ny),tau(1:nx,1:ny)
   real*8,intent(out) :: enin(1:nx,1:ny),Ener(1:nx,1:ny)
   integer(KNDI),intent(out) :: EOS(1:nx,1:ny)
   real*8,intent(out) :: vertx(0:nx,0:ny),verty(0:nx,0:ny)
   character(lenofchar),intent(out) :: node_type_char(0:nx,0:ny)
   real(kndr),intent(out)           :: node_type_direction(0:nx,0:ny,2)
   real(kndr),intent(out)           :: node_type_velocity(0:nx,0:ny,2)
   real(kndr),intent(out)           :: node_type_pressure(0:nx,0:ny) 

   integer :: ix, iy,i,j
   real*8 :: r_left,r_rght,theta_left,theta_rght
   real*8 :: dr,dtheta,r,theta,tmp_theta(0:6)
   
   print*, '---> Initial_Mesh_Phys_Sod_Sector_NoEqual_angle_0'
   
   r_left = 0.d0
   r_rght = 1.d0
   theta_left = 0.d0
   theta_rght = PI/2.d0
   
   dr = (r_rght - r_left)/nx
   dtheta = (theta_rght-theta_left)/ny

   tmp_theta(0)  = 0.d0/180.d0*PI   
   tmp_theta(1)  = 4.1d0/180.d0*PI
   tmp_theta(2)  = 10.3d0/180.d0*PI
   tmp_theta(3)  = 33.1d0/180.d0*PI
   tmp_theta(4)  = 51.6d0/180.d0*PI
   tmp_theta(5)  = 82.4d0/180.d0*PI
   tmp_theta(ny) = 90.d0/180.d0*PI

   do iy = 0 , ny
      theta = tmp_theta(iy)
      do ix = 0 , nx
         r = r_left + ix * dr
         vertx(ix,iy) = r * cos(theta)
         verty(ix,iy) = r * sin(theta)
         if(iy==0) then        !> x axis
            vertx(ix,iy) = r 
            verty(ix,iy) = 0.d0
         elseif(iy==ny) then   !> y axis
            vertx(ix,iy) = 0.d0 
            verty(ix,iy) = r
         endif
      enddo
   enddo
   
   !> Phys
   do  ix = 1 , nx
   do  iy = 1 , ny
      if(ix <= nx/2) then
         den(ix,iy)  = 1.d0
         velx(ix,iy) = 0.d0
         vely(ix,iy) = 0.d0
         pre(ix,iy)  = 1.d0
         EOS(ix,iy)  = eos_num
      else
         den(ix,iy)  = 0.125d0
         velx(ix,iy) = 0.d0
         vely(ix,iy) = 0.d0
         pre(ix,iy)  = 0.1d0
         EOS(ix,iy)  = eos_num
      endif
   enddo
   enddo
   
   do  ix = 1 , nx
   do  iy = 1 , ny
      if(ix <= nx/2) then
         Call General_EOS(EOS(ix,iy))
         sound(ix,iy) = sqrt( gammar*pre(ix,iy)/den(ix,iy) )
         enin(ix,iy)  = pre(ix,iy)/(gammar-1.d0)/den(ix,iy)
         tau(ix,iy)   = 1.d0 / den(ix,iy)
         Ener(ix,iy)  = enin(ix,iy) + &
                        ((velx(ix,iy))**2+(vely(ix,iy))**2)/2.d0
      else
         Call General_EOS(EOS(ix,iy))
         sound(ix,iy) = sqrt( gammar*pre(ix,iy)/den(ix,iy) )
         enin(ix,iy)  = pre(ix,iy)/(gammar-1.d0)/den(ix,iy)
         tau(ix,iy)   = 1.d0 / den(ix,iy)
         Ener(ix,iy)  = enin(ix,iy) + &
                        ((velx(ix,iy))**2+(vely(ix,iy))**2)/2.d0
      endif
   enddo
   enddo 
   
   !> bndy condition 
   !> x_bgn
   i = 0
   do j = 1 , ny - 1
      node_type_char(i,j)    = 'fix'
      node_type_direction(i,j,1) = 0.d0      
      node_type_direction(i,j,2) = 0.d0      
      node_type_velocity(i,j,1)  = 0.d0
      node_type_velocity(i,j,2)  = 0.d0
      node_type_pressure(i,j)     = 0.d0
   enddo 
   !> x_end
   i = nx
   do j = 1 , ny - 1
      node_type_char(i,j)    = 'fix'
      node_type_direction(i,j,1) = 0.d0      
      node_type_direction(i,j,2) = 0.d0      
      node_type_velocity(i,j,1)  = 0.d0
      node_type_velocity(i,j,2)  = 0.d0
      node_type_pressure(i,j)     = 0.d0
   enddo    
   !> y_bgn
   j = 0
   do i = 1 , nx - 1
      node_type_char(i,j)    = 'slide_line'
      node_type_direction(i,j,1) = 1.d0      
      node_type_direction(i,j,2) = 0.d0      
      node_type_velocity(i,j,1)  = 0.d0
      node_type_velocity(i,j,2)  = 0.d0
      node_type_pressure(i,j)     = 0.d0
   enddo 
   !> y_end
   j = ny
   do i = 1 , nx - 1
      node_type_char(i,j)    = 'slide_line'
      node_type_direction(i,j,1) = 0.d0      
      node_type_direction(i,j,2) = 1.d0      
      node_type_velocity(i,j,1)  = 0.d0
      node_type_velocity(i,j,2)  = 0.d0
      node_type_pressure(i,j)     = 0.d0
   enddo    

   !> x_bgn,y_bgn
   i = 0 ; j = 0
   node_type_char(i,j)    = 'fix'
   node_type_direction(i,j,1) = 0.d0      
   node_type_direction(i,j,2) = 0.d0      
   node_type_velocity(i,j,1)  = 0.d0
   node_type_velocity(i,j,2)  = 0.d0
   node_type_pressure(i,j)     = 0.d0   
   !> x_bgn,y_end
   i = 0 ; j = ny
   node_type_char(i,j)    = 'fix'
   node_type_direction(i,j,1) = 0.d0      
   node_type_direction(i,j,2) = 0.d0      
   node_type_velocity(i,j,1)  = 0.d0
   node_type_velocity(i,j,2)  = 0.d0
   node_type_pressure(i,j)     = 0.d0   
   !> x_end,y_bgn
   i = nx ; j = 0
   node_type_char(i,j)    = 'fix'
   node_type_direction(i,j,1) = 0.d0      
   node_type_direction(i,j,2) = 0.d0      
   node_type_velocity(i,j,1)  = 0.d0
   node_type_velocity(i,j,2)  = 0.d0
   node_type_pressure(i,j)     = 0.d0   
   !> x_end,y_end
   i = nx ; j = ny
   node_type_char(i,j)    = 'fix'
   node_type_direction(i,j,1) = 0.d0      
   node_type_direction(i,j,2) = 0.d0      
   node_type_velocity(i,j,1)  = 0.d0
   node_type_velocity(i,j,2)  = 0.d0
   node_type_pressure(i,j)     = 0.d0   
      

   print*, '<--- Initial_Mesh_Phys_Sod_Sector_NoEqual_angle_0'
   
end subroutine Initial_Mesh_Phys_Sod_Sector_NoEqual_angle_0

End module mdu_sod_2D_cat_cyl