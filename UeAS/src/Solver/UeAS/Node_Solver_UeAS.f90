
Subroutine Node_Solver_Cat_UeAS(Ini_CDT,MeshCY2D,PhysCY2D)
   use COMMON_TYPE_Condition
   use COMMON_TYPE_Mesh
   use COMMON_TYPE_Phys
   use COMMON_XXX_Constant
   implicit none

   Type(COMMON_TYPE_CDT),intent(inout) :: Ini_CDT
   Type(Mesh_Cyl_2D),intent(inout) :: MeshCY2D
   Type(Phys_Cyl_2D),intent(inout) :: PhysCY2D

   Call Node_Solver_Cat_UeAS_InnerNode_add_pian(Ini_CDT%nx,Ini_CDT%ny,PhysCY2D%den,PhysCY2D%vel_x, &
         PhysCY2D%vel_y,PhysCY2D%pre,PhysCY2D%sound,PhysCY2D%Vstar_x,PhysCY2D%Vstar_y,     &
         PhysCY2D%PstarC1,MeshCY2D%vertex_x,MeshCY2D%vertex_y)

   Call Node_Solver_Cat_UeAS_BndyNode(Ini_CDT,MeshCY2D,PhysCY2D)         

End Subroutine Node_Solver_Cat_UeAS


Subroutine Node_Solver_Cat_UeAS_InnerNode_add_pian(nx,ny,den,velx,vely,pre,sound,Vstarx,Vstary,Pstar1,&
            vertx,verty)
   use COMMON_XXX_Constant
   implicit none
   
   integer,intent(in) :: nx,ny
   Real*8,intent(in)  :: den(1:nx,1:ny)
   Real*8,intent(in)  :: pre(1:nx,1:ny),sound(1:nx,1:ny)
   Real*8,intent(in)  :: velx(1:nx,1:ny),vely(1:nx,1:ny)
   Real*8,intent(in)  :: vertx(0:nx,0:ny),verty(0:nx,0:ny)
   
   Real*8,intent(out) :: Vstarx(0:nx,0:ny),Vstary(0:nx,0:ny)
   Real*8,intent(out) :: Pstar1(0:nx,0:ny)
    
   integer :: i,j,i2,j2,i3,j3,i4,j4,side
   real*8  :: temptA,temptb,temptc,temptd,smx,smy,starvk,vv,phi
   real*8  :: alength,anormalx,anormaly,anormalxnn,anormalynn
   real*8  :: x1,y1,x2,y2,x3,y3,addalpha,rx2,rx3,ry2,ry3,phi2,phi3,delphi
   real*8  :: njr2,njr3,Ajr11,Ajr12,Ajr22,right_hand,Gama_r,pre_r
   real*8  :: Vjr11,Vjr12,Vjr21,Vjr22,vx,vy,aa,bb,cc,theta_m,theta_d,u_bar
   real*8  :: rho1,sod1,pre1,vel1,ustar(4),pstar(4),cos2theta,costheta,sintheta
   real*8  :: rho2,sod2,pre2,vel2,ustar_max,pstar_max,l6,l7,alpha4,vr1,vr2
   real*8  :: area,l2,l3,l4,l5
   
   do j=1,ny-1
   do i=1,nx-1

      temptA=0.d0
      temptB=0.d0
      temptC=0.d0
      temptD=0.d0
      SMx=0.d0
      SMy=0.d0
      right_hand = 0.d0
      Gama_r = 0.d0         
         
      do side = 1 , 4      !> 4 means quadrangle
         select case(side)
         case(1)
            i2 = i+1 ; j2 = j
            i3 = i ; j3 = j + 1
            i4 = i+1 ; j4 = j+1
         case(2)
            i2 = i ; j2 = j+1
            i3 = i-1 ; j3 = j
            i4 = i ; j4 = j+1
         case(3)
            i2 = i-1 ; j2 = j
            i3 = i ; j3 = j-1
            i4 = i ; j4 = j     
         case(4)
            i2 = i ; j2 = j-1
            i3 = i+1 ; j3 = j
            i4 = i+1 ; j4 = j   
         case default
            stop 'side is error, pls call xihua.'
         endselect      
   
         x1=vertx(i,j)   ; y1=verty(i,j)
         x2=vertx(i2,j2) ; y2=verty(i2,j2)
         rx2 = x2 - x1    ; ry2 = y2 - y1
         Call Cal_theta(rx2,ry2,phi2)

         x3=vertx(i3,j3) ; y3=verty(i3,j3)
         rx3 = x3 - x1    ; ry3 = y3 - y1
         Call Cal_theta(rx3,ry3,phi3)
         
         Call Cal_njr(phi2,phi3,njr2,njr3)
         if(phi3<phi2) then
            delphi = phi3+2.d0*PI-phi2
         else
            delphi = phi3 - phi2
         endif    
                
         Ajr11 = (2.d0*delphi+dsin(2.d0*phi3)-dsin(2.d0*phi2))*den(i4,j4)*sound(i4,j4) / 4.d0
         Ajr22 = (2.d0*delphi-dsin(2.d0*phi3)+dsin(2.d0*phi2))*den(i4,j4)*sound(i4,j4) / 4.d0
         Ajr12 = (dcos(2.d0*phi2)-dcos(2.d0*phi3))*den(i4,j4)*sound(i4,j4) / 4.d0
         temptA = temptA + Ajr11
         temptB = temptB + Ajr22
         temptC = temptC + Ajr12         
         
         phi = phi2 + delphi/2.d0
         Vjr11 = cos(phi)*cos(phi)*den(i4,j4)*sound(i4,j4)
         Vjr22 = sin(phi)*sin(phi)*den(i4,j4)*sound(i4,j4)
         Vjr12 = cos(phi)*sin(phi)*den(i4,j4)*sound(i4,j4)
               
         SMx = SMx + Vjr11*velx(i4,j4) + Vjr12*vely(i4,j4) - pre(i4,j4)*cos(phi)*sin(delphi/2.d0)
         SMy = SMy + Vjr12*velx(i4,j4) + Vjr22*vely(i4,j4) - pre(i4,j4)*sin(phi)*sin(delphi/2.d0) 
         
         phi = phi2 + delphi/2.d0
         Gama_r = Gama_r + 1.d0/sound(i4,j4)/den(i4,j4)
         right_hand = right_hand + 1.d0/sound(i4,j4)/den(i4,j4)*pre(i4,j4) - &
                     (velx(i4,j4)*(cos(phi)) + vely(i4,j4)*(sin(phi))) / sin(delphi/2.d0)           

      enddo
      Vstarx(i,j)=(SMx*temptB-SMy*temptC)/(temptA*temptB-temptC*temptC) 
      Vstary(i,j)=(SMy*temptA-SMx*temptC)/(temptA*temptB-temptC*temptC)      
      Pstar1(i,j) = right_hand / Gama_r

      area = 0.d0    
      do side = 1 , 4      !> 4 means quadrangle
         select case(side)
         case(1)
            i2 = i+1 ; j2 = j
            i3 = i ; j3 = j + 1
            i4 = i+1 ; j4 = j+1            
         case(2)
            i2 = i ; j2 = j+1
            i3 = i-1 ; j3 = j
            i4 = i ; j4 = j+1
         case(3)
            i2 = i-1 ; j2 = j
            i3 = i ; j3 = j-1
            i4 = i ; j4 = j     
         case(4)
            i2 = i ; j2 = j-1
            i3 = i+1 ; j3 = j
            i4 = i+1 ; j4 = j   
         case default         
            stop 'side is error, pls call xihua.'
         endselect      
   
         x1=vertx(i,j)   ; y1=verty(i,j)
         x2=vertx(i2,j2) ; y2=verty(i2,j2)
         rx2 = x2 - x1   ; ry2 = y2 - y1
         Call Cal_theta(rx2,ry2,phi2)

         x3=vertx(i3,j3) ; y3=verty(i3,j3)
         rx3 = x3 - x1   ; ry3 = y3 - y1
         Call Cal_theta(rx3,ry3,phi3)
         
         l2 = dsqrt(rx2**2 + ry2**2)
         l3 = dsqrt(rx3**2 + ry3**2)
         
         Call Cal_njr(phi2,phi3,njr2,njr3)
         delphi = phi3 - phi2   
         if(delphi < 0.d0) then
            phi2 = phi2 - 2.d0*PI
            delphi = phi3 - phi2
         endif         
         area = area + abs(l2*l3*dsin(delphi))               
      enddo
      
      l4 = dsqrt ( ( vertx(i+1,j) - vertx(i-1,j) )**2 + &
                  ( verty(i+1,j) - verty(i-1,j) )**2 )
      l5 = dsqrt ( ( vertx(i,j+1) - vertx(i,j-1) )**2 + &
                  ( verty(i,j+1) - verty(i,j-1) )**2 )
      
      cos2theta = (area/(l4*l5))
      
      l6 = dsqrt ( ( vertx(i,j) - vertx(i,j+1) )**2 + &
                  ( verty(i,j) - verty(i,j+1) )**2 )
      l7 = dsqrt ( ( vertx(i,j) - vertx(i,j-1) )**2 + &
                  ( verty(i,j) - verty(i,j-1) )**2 )      
      if(l6 > l7) then
         costheta = dsqrt( (1.d0 + cos2theta)/2.d0 )
         sintheta = dsqrt( (1.d0 - cos2theta)/2.d0 )
      else
         costheta =  dsqrt( (1.d0 + cos2theta)/2.d0 )
         sintheta = -dsqrt( (1.d0 - cos2theta)/2.d0 )          
      endif
      if(abs(sintheta) < eps) then
         sintheta = 0.d0
         costheta = 1.d0
      endif  
      vx = Vstarx(i,j) * costheta   + Vstary(i,j) * sintheta
      vy = Vstarx(i,j) *(-sintheta) + Vstary(i,j) * costheta 
      if( sqrt(vx*vx+vy*vy) < 1.d-10) then
         Vstarx(i,j) = 0.d0
         Vstary(i,j) = 0.d0 
         cycle
      endif      
      
      theta_m = acos(vx/sqrt(vx*vx+vy*vy))

      right_hand = 0.d0
      Gama_r     = 0.d0  
      do side = 1 , 4      !> 4 means quadrangle
         select case(side)
         case(1)
            i2 = i+1 ; j2 = j
            i3 = i ; j3 = j + 1
            i4 = i+1 ; j4 = j+1 
         case(2)
            i2 = i ; j2 = j+1
            i3 = i-1 ; j3 = j
            i4 = i ; j4 = j+1
         case(3)
            i2 = i-1 ; j2 = j
            i3 = i ; j3 = j-1
            i4 = i ; j4 = j   
         case(4)
            i2 = i ; j2 = j-1
            i3 = i+1 ; j3 = j
            i4 = i+1 ; j4 = j   
         case default
            stop 'side is error, pls call xihua.'
         endselect
         
         x1=vertx(i,j)   ; y1=verty(i,j)
         x2=vertx(i2,j2) ; y2=verty(i2,j2)
         rx2 = x2 - x1    ; ry2 = y2 - y1
         Call Cal_theta(rx2,ry2,phi2)

         x3=vertx(i3,j3) ; y3=verty(i3,j3)
         rx3 = x3 - x1    ; ry3 = y3 - y1
         Call Cal_theta(rx3,ry3,phi3)
         
         Call Cal_njr(phi2,phi3,njr2,njr3)
         if(phi3<phi2) then
            delphi = phi3+2.d0*PI-phi2
         else
            delphi = phi3 - phi2
         endif
         
         phi = phi2 + delphi/2.d0 
         Gama_r = Gama_r + (vx / sqrt(vx*vx+vy*vy)*cos(phi) + vy / sqrt(vx*vx+vy*vy)*sin(phi))
         right_hand = right_hand + (pre(i4,j4)-Pstar1(i,j))/sound(i4,j4)/den(i4,j4)*sin(delphi/2.d0) - &
                     (velx(i4,j4)*(cos(phi)) + vely(i4,j4)*(sin(phi)))
      enddo            
      u_bar = right_hand / Gama_r      
      Vstarx(i,j) = abs(u_bar) * vx / sqrt(vx*vx+vy*vy)
      Vstary(i,j) = abs(u_bar) * vy / sqrt(vx*vx+vy*vy)
      
	enddo
	enddo
      
end subroutine Node_Solver_Cat_UeAS_InnerNode_add_pian

