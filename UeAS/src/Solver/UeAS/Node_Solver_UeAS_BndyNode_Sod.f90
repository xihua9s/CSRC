
Subroutine Node_Solver_Cat_UeAS_BndyNode_Sod_Sector(nx,ny,den,velx,vely,pre,sound, &
            Vstarx,Vstary,Pstar1,vertx,verty,node_type_char)
   use COMMON_TYPE_Condition
   use COMMON_XXX_Constant
   implicit none
   
   integer,intent(in) :: nx,ny
   Real*8,intent(in)  :: den(1:nx,1:ny)
   Real*8,intent(in)  :: pre(1:nx,1:ny),sound(1:nx,1:ny)
   Real*8,intent(in)  :: velx(1:nx,1:ny),vely(1:nx,1:ny)
   Real*8,intent(in)  :: vertx(0:nx,0:ny),verty(0:nx,0:ny)
   character(lenofchar),intent(in) :: node_type_char(0:nx,0:ny)
   
   Real*8,intent(out) :: Vstarx(0:nx,0:ny),Vstary(0:nx,0:ny)
   Real*8,intent(out) :: Pstar1(0:nx,0:ny)

   integer :: i,j,i2,j2,i3,j3,i4,j4,side
   real*8  :: temptA,temptb,temptc,smx,smy,starvk,vv,temptAhat,temptBhat,temptChat
   real*8  :: alength,anormalx,anormaly,anormalxnn,anormalynn
   real*8  :: x1,y1,x2,y2,x3,y3,addalpha,rx2,rx3,ry2,ry3,phi2,phi3,delphi
   real*8  :: njr2,njr3,Ajr11,Ajr12,Ajr22,right_hand,Gama_r,pre_r
   real*8  :: normal_wall_x,normal_wall_y,SMxhat,SMyhat,tmpl1,tmpc1,tmpv
   real*8  :: l2,l3,l4,l5,tmp_theta,tmp_theta_1,area,l22,l33,lc,sigma,theta
   real*8  :: rho1,sod1,pre1,vel1,ustar(4),pstar(4)
   real*8  :: rho2,sod2,pre2,vel2,ustar_max,pstar_max   

   !> Bottom
   j = 0
   do i=1,nx-1
      
      i2 = i ; j2 = j+1
      i3 = i+1 ; j3 = j+1
      i4 = i ; j4 = j+1
         
      x1=vertx(i,j)   ; y1=verty(i,j)
      x2=vertx(i2,j2) ; y2=verty(i2,j2)         
      Call Normal(x1,y1,x2,y2,anormalx,anormaly)
      
      rho1 = den(i3,j3)
      sod1 = sound(i3,j3)
      vel1 = velx(i3,j3)*anormalx + vely(i3,j3)*anormaly
      pre1 = pre(i3,j3)
      
      rho2 = den(i4,j4)
      sod2 = sound(i4,j4)
      vel2 = velx(i4,j4)*anormalx + vely(i4,j4)*anormaly
      pre2 = pre(i4,j4)
      
      ustar_max = (rho1*sod1*vel1 + rho2*sod2*vel2)/(rho1*sod1 + rho2*sod2) + &
                    (pre1 - pre2)/(rho1*sod1 + rho2*sod2)
      pstar_max = (pre1*rho2*sod2 + pre2*rho1*sod1)/(rho1*sod1 + rho2*sod2) + &
                  (vel1-vel2)*(rho1*sod1*rho2*sod2)/(rho1*sod1 + rho2*sod2) 
      
      Pstar1(i,j) = pstar_max
      Vstarx(i,j) = -ustar_max
      Vstary(i,j) = 0.d0      

	enddo

   !> top 
   j = ny
   do i=1,nx-1     
      i2 = i ; j2 = j-1
      i3 = i ; j3 = j
      i4 = i+1 ; j4 = j   
         
      x1=vertx(i,j)   ; y1=verty(i,j)
      x2=vertx(i2,j2) ; y2=verty(i2,j2)         
      Call Normal(x1,y1,x2,y2,anormalx,anormaly)
      
      rho1 = den(i3,j3)
      sod1 = sound(i3,j3)
      vel1 = velx(i3,j3)*anormalx + vely(i3,j3)*anormaly
      pre1 = pre(i3,j3)
      
      rho2 = den(i4,j4)
      sod2 = sound(i4,j4)
      vel2 = velx(i4,j4)*anormalx + vely(i4,j4)*anormaly
      pre2 = pre(i4,j4)      
   
      ustar_max = (rho1*sod1*vel1 + rho2*sod2*vel2)/(rho1*sod1 + rho2*sod2) + &
                    (pre1 - pre2)/(rho1*sod1 + rho2*sod2)

      pstar_max = (pre1*rho2*sod2 + pre2*rho1*sod1)/(rho1*sod1 + rho2*sod2) + &
                  (vel1-vel2)*(rho1*sod1*rho2*sod2)/(rho1*sod1 + rho2*sod2) 
      
      Pstar1(i,j) = pstar_max
      Vstary(i,j) = ustar_max      
      Vstarx(i,j) = 0.d0
     
   enddo
   
   !> Left Free Boundary (i==0)
   i = 0
   do j = 1 , ny-1
      Vstarx(i,j) = 0.d0 
      Vstary(i,j) = 0.d0 
      Pstar1(i,j) = pre(i+1,j)
   enddo

   !> Right Free Boundary (i==nx)
   i = nx
   do j = 1 , ny-1
      Vstarx(i,j) = 0.d0
      Vstary(i,j) = 0.d0
      Pstar1(i,j) = pre(i,j)
   enddo

   !> Left + Bottom Node (i==0,j==0)
   i = 0 ; j = 0
   Vstarx(i,j) = 0.d0 
   Vstary(i,j) = 0.d0 
   Pstar1(i,j) = pre(i+1,j+1)
      
   !> Left + Top Node (i==0,j==ny)
   i = 0 ; j = ny
   Vstarx(i,j) = 0.d0 
   Vstary(i,j) = 0.d0 
   Pstar1(i,j) = pre(i+1,j)

   !> Right + Bottom Node (i==nx,j==0)
   i = nx ; j = 0
   Vstarx(i,j) = 0.d0 
   Vstary(i,j) = 0.d0 
   Pstar1(i,j) = pre(i,j+1)

   !> Right + Top Node (i==nx,j==ny)
   i = nx ; j = ny
   Vstarx(i,j) = 0.d0 
   Vstary(i,j) = 0.d0 
   Pstar1(i,j) = pre(i,j)     
   
End Subroutine Node_Solver_Cat_UeAS_BndyNode_Sod_Sector


