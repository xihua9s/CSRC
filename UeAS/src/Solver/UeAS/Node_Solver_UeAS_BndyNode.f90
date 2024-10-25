

Subroutine Node_Solver_Cat_UeAS_BndyNode(Ini_CDT,MeshCY2D,PhysCY2D)
   use COMMON_TYPE_Condition
   use COMMON_TYPE_Mesh
   use COMMON_TYPE_Phys
   use COMMON_XXX_Constant
   implicit none

   Type(COMMON_TYPE_CDT),intent(inout) :: Ini_CDT
   Type(Mesh_Cyl_2D),intent(inout) :: MeshCY2D
   Type(Phys_Cyl_2D),intent(inout) :: PhysCY2D
   
   Call Node_Solver_Cat_UeAS_BndyNode_add_pian(Ini_CDT,Ini_CDT%nx,Ini_CDT%ny,PhysCY2D%den,PhysCY2D%vel_x, &
         PhysCY2D%vel_y,PhysCY2D%pre,PhysCY2D%sound,PhysCY2D%Vstar_x,PhysCY2D%Vstar_y,    &
         PhysCY2D%PstarC1,MeshCY2D%vertex_x,MeshCY2D%vertex_y,MeshCY2D%node_type_char,&
         MeshCY2D%node_type_direction,MeshCY2D%node_type_velocity,MeshCY2D%node_type_pressure,&
         Ini_CDT%Gou_Xing,Ini_CDT%Origin_YesNo)  
   
End Subroutine Node_Solver_Cat_UeAS_BndyNode  


Subroutine Node_Solver_Cat_UeAS_BndyNode_add_pian(Ini_CDT,nx,ny,den,velx,vely,pre,&
            sound,Vstarx,Vstary,Pstar1,vertx,verty,node_type_char,     &
            node_type_direction,node_type_velocity,node_type_pressure, &
            Gou_Xing,Origin_YesNo)
   use COMMON_TYPE_Condition
   use COMMON_XXX_Constant
   implicit none
   
   Type(COMMON_TYPE_CDT),intent(inout) :: Ini_CDT   
   integer,intent(in) :: nx,ny
   Real*8,intent(in)  :: den(1:nx,1:ny)
   Real*8,intent(in)  :: pre(1:nx,1:ny),sound(1:nx,1:ny)
   Real*8,intent(in)  :: velx(1:nx,1:ny),vely(1:nx,1:ny)
   Real*8,intent(in)  :: vertx(0:nx,0:ny),verty(0:nx,0:ny)
   character(lenofchar),intent(in) :: node_type_char(0:nx,0:ny)
   real(kndr),intent(in)           :: node_type_direction(0:nx,0:ny,2)
   real(kndr),intent(in)           :: node_type_velocity(0:nx,0:ny,2)
   real(kndr),intent(in)           :: node_type_pressure(0:nx,0:ny)
   character(36),intent(in)        :: Gou_Xing,Origin_YesNo
   
   Real*8,intent(inout) :: Vstarx(0:nx,0:ny),Vstary(0:nx,0:ny)
   Real*8,intent(inout) :: Pstar1(0:nx,0:ny)
   
   integer :: i,j,i2,j2,i3,j3,i4,j4,side
   real*8  :: temptA,temptb,temptc,temptd,smx,smy,starvk,vv,phi
   real*8  :: alength,anormalx,anormaly,anormalxnn,anormalynn
   real*8  :: x1,y1,x2,y2,x3,y3,addalpha,rx2,rx3,ry2,ry3,phi2,phi3,delphi
   real*8  :: njr2,njr3,Ajr11,Ajr12,Ajr22,right_hand,Gama_r,pre_r
   real*8  :: Vjr11,Vjr12,Vjr22,vx,vy,aa,bb,cc,theta_m,theta_d,u_bar
   real*8  :: rho1,sod1,pre1,vel1,ustar(4),pstar(4),cos2theta,costheta,sintheta
   real*8  :: rho2,sod2,pre2,vel2,ustar_max,pstar_max,l6,l7,alpha4,vr1,vr2
   real*8  :: area,l2,l3,l4,l5,Ione(2,2),tangent(2),tot(2,2),tot_A_tot(2,2),barA(2,2),barB(2)
   real(kndr)      :: norm(2),tvel(2),tpre,A(2,2),B(2)
      
   Ione(1,1) = 1.d0 
   Ione(1,2) = 0.d0 
   Ione(2,1) = 0.d0 
   Ione(2,2) = 1.d0    

   !> x_bgn
   i = 0
   Do j = 1 , ny - 1
      tangent = node_type_direction(i,j,1:2)
      tvel    = node_type_velocity(i,j,1:2)
      tpre    = node_type_pressure(i,j)    
      
      if(Gou_Xing == 'Quadrangle' .or. Origin_YesNo == 'No') then

         Select case(node_type_char(i,j))
         
         Case ('fix')
            Vstarx(i,j) = 0.d0  
            Vstary(i,j) = 0.d0 
            
            right_hand  = 0.d0
            Gama_r      = 0.d0         
            do side = 1 , 4      !> 4 means quadrangle
               select case(side)
               case(1)
                  i2 = i+1 ; j2 = j
                  i3 = i ; j3 = j + 1
                  i4 = i+1 ; j4 = j+1
               case(2)
                  ! i2 = i ; j2 = j+1
                  ! i3 = i-1 ; j3 = j
                  ! i4 = i ; j4 = j+1
                  cycle
               case(3)
                  ! i2 = i-1 ; j2 = j
                  ! i3 = i ; j3 = j-1
                  ! i4 = i ; j4 = j     
                  cycle
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
               Gama_r = Gama_r + 1.d0/sound(i4,j4)/den(i4,j4)
               right_hand = right_hand + 1.d0/sound(i4,j4)/den(i4,j4)*pre(i4,j4) - &
                           (velx(i4,j4)*(cos(phi)) + vely(i4,j4)*(sin(phi))) / sin(delphi/2.d0)
            enddo       
            Pstar1(i,j) = right_hand / Gama_r     

         Case ('given_velocity')      
            Vstarx(i,j) = tvel(1)  
            Vstary(i,j) = tvel(2)
            
            right_hand  = 0.d0
            Gama_r      = 0.d0         
            do side = 1 , 4      !> 4 means quadrangle
               select case(side)
               case(1)
                  i2 = i+1 ; j2 = j
                  i3 = i ; j3 = j + 1
                  i4 = i+1 ; j4 = j+1
               case(2)
                  ! i2 = i ; j2 = j+1
                  ! i3 = i-1 ; j3 = j
                  ! i4 = i ; j4 = j+1
                  cycle
               case(3)
                  ! i2 = i-1 ; j2 = j
                  ! i3 = i ; j3 = j-1
                  ! i4 = i ; j4 = j     
                  cycle
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
               Gama_r = Gama_r + 1.d0/sound(i4,j4)/den(i4,j4)
               right_hand = right_hand + 1.d0/sound(i4,j4)/den(i4,j4)*pre(i4,j4) - &
                           (velx(i4,j4)*(cos(phi)) + vely(i4,j4)*(sin(phi))) / sin(delphi/2.d0)
            enddo       
            Pstar1(i,j) = right_hand / Gama_r       
         
         Case ('given_pressure')         
            stop 'Not Ready, 189xihua'            
         Case ('slide_line')
         
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
                  ! i2 = i ; j2 = j+1
                  ! i3 = i-1 ; j3 = j
                  ! i4 = i ; j4 = j+1
                  cycle
               case(3)
                  ! i2 = i-1 ; j2 = j
                  ! i3 = i ; j3 = j-1
                  ! i4 = i ; j4 = j     
                  cycle
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
                      
               !> 2022-06-16
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

            A(1,1) = temptA
            A(1,2) = temptC
            A(2,1) = temptC
            A(2,2) = temptB
            B(1)   = SMx
            B(2)   = SMy
            call tangent_otimes_tangent(tangent,tot)
            call tot_dot_A_dot_tot(tot,A,tot_A_tot)
            barA = tot_A_tot + (temptA+temptB)*(Ione - tot)            
            call tot_dot_B(tot,B,barB)
            vx = (barB(1)*barA(2,2)-barB(2)*barA(1,2))/(barA(1,1)*barA(2,2)-barA(1,2)*barA(1,2))
            vy = (barB(2)*barA(1,1)-barB(1)*barA(1,2))/(barA(1,1)*barA(2,2)-barA(1,2)*barA(1,2))
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
                  ! i2 = i ; j2 = j+1
                  ! i3 = i-1 ; j3 = j
                  ! i4 = i ; j4 = j+1
                  cycle
               case(3)
                  ! i2 = i-1 ; j2 = j
                  ! i3 = i ; j3 = j-1
                  ! i4 = i ; j4 = j     
                  cycle
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

         Case default
            
            print*, 'node_type_char(i,j)=',node_type_char(i,j)
            stop 'node_type_char(i,j) is error, pls call xihua f8a'
      
         End Select

      elseif(Gou_Xing == 'Sector' .and. Origin_YesNo == 'Yes') then

         Select case(node_type_char(i,j))
         
         Case ('fix')
            Vstarx(i,j) = 0.d0  
            Vstary(i,j) = 0.d0 
            
            right_hand  = 0.d0
            Gama_r      = 0.d0         
            do side = 1 , 4      !> 4 means quadrangle
               select case(side)
               case(1)
                  i2 = i+1 ; j2 = j
                  i3 = i+1 ; j3 = j + 1
                  i4 = i+1 ; j4 = j+1
               case(2)
                  ! i2 = i ; j2 = j+1
                  ! i3 = i-1 ; j3 = j
                  ! i4 = i ; j4 = j+1
                  cycle
               case(3)
                  ! i2 = i-1 ; j2 = j
                  ! i3 = i ; j3 = j-1
                  ! i4 = i ; j4 = j     
                  cycle
               case(4)
                  i2 = i+1 ; j2 = j-1
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
               Gama_r = Gama_r + 1.d0/sound(i4,j4)/den(i4,j4)
               right_hand = right_hand + 1.d0/sound(i4,j4)/den(i4,j4)*pre(i4,j4) - &
                           (velx(i4,j4)*(cos(phi)) + vely(i4,j4)*(sin(phi)))
            enddo       
            Pstar1(i,j) = right_hand / Gama_r     
         Case ('given_velocity')      
            Vstarx(i,j) = tvel(1)  
            Vstary(i,j) = tvel(2)
            
            right_hand  = 0.d0
            Gama_r      = 0.d0         
            do side = 1 , 4      !> 4 means quadrangle
               select case(side)
               case(1)
                  i2 = i+1 ; j2 = j
                  i3 = i+1 ; j3 = j + 1
                  i4 = i+1 ; j4 = j+1
               case(2)
                  ! i2 = i ; j2 = j+1
                  ! i3 = i-1 ; j3 = j
                  ! i4 = i ; j4 = j+1
                  cycle
               case(3)
                  ! i2 = i-1 ; j2 = j
                  ! i3 = i ; j3 = j-1
                  ! i4 = i ; j4 = j     
                  cycle
               case(4)
                  i2 = i+1 ; j2 = j-1
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
               Gama_r = Gama_r + 1.d0/sound(i4,j4)/den(i4,j4)
               right_hand = right_hand + 1.d0/sound(i4,j4)/den(i4,j4)*pre(i4,j4) - &
                           (velx(i4,j4)*(cos(phi)) + vely(i4,j4)*(sin(phi))) / sin(delphi/2.d0)
            enddo       
            Pstar1(i,j) = right_hand / Gama_r       
         
         Case ('given_pressure')
            stop 'Not Ready, 459xihua'                     
         Case ('slide_line')
         
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
                  i3 = i+1 ; j3 = j + 1
                  i4 = i+1 ; j4 = j+1
               case(2)
                  ! i2 = i ; j2 = j+1
                  ! i3 = i-1 ; j3 = j
                  ! i4 = i ; j4 = j+1
                  cycle
               case(3)
                  ! i2 = i-1 ; j2 = j
                  ! i3 = i ; j3 = j-1
                  ! i4 = i ; j4 = j     
                  cycle
               case(4)
                  i2 = i+1 ; j2 = j-1
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
                      
               !> 2022-06-16
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

            A(1,1) = temptA
            A(1,2) = temptC
            A(2,1) = temptC
            A(2,2) = temptB
            B(1)   = SMx
            B(2)   = SMy
            call tangent_otimes_tangent(tangent,tot)
            call tot_dot_A_dot_tot(tot,A,tot_A_tot)
            barA = tot_A_tot + (temptA+temptB)*(Ione - tot)            
            call tot_dot_B(tot,B,barB)
            vx = (barB(1)*barA(2,2)-barB(2)*barA(1,2))/(barA(1,1)*barA(2,2)-barA(1,2)*barA(1,2))
            vy = (barB(2)*barA(1,1)-barB(1)*barA(1,2))/(barA(1,1)*barA(2,2)-barA(1,2)*barA(1,2))
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
                  i3 = i+1 ; j3 = j + 1
                  i4 = i+1 ; j4 = j+1
               case(2)
                  ! i2 = i ; j2 = j+1
                  ! i3 = i-1 ; j3 = j
                  ! i4 = i ; j4 = j+1
                  cycle
               case(3)
                  ! i2 = i-1 ; j2 = j
                  ! i3 = i ; j3 = j-1
                  ! i4 = i ; j4 = j     
                  cycle
               case(4)
                  i2 = i+1 ; j2 = j-1
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
         Case default
            
            print*, 'node_type_char(i,j)=',node_type_char(i,j)
            stop 'node_type_char(i,j) is error, pls call xihua f8a'
      
         End Select

      else 
         stop 'xihuafd9a'         
      endif 
   Enddo
    
   !> x_end
   i = nx
   Do j = 1 , ny - 1

      tangent = node_type_direction(i,j,1:2)
      tvel    = node_type_velocity(i,j,1:2)
      tpre    = node_type_pressure(i,j)    
      Select case(node_type_char(i,j))
      
      Case ('fix')
         Vstarx(i,j) = 0.d0  
         Vstary(i,j) = 0.d0 
         
         right_hand  = 0.d0
         Gama_r      = 0.d0         
         do side = 1 , 4      !> 4 means quadrangle
            select case(side)
            case(1)
               ! i2 = i+1 ; j2 = j
               ! i3 = i ; j3 = j + 1
               ! i4 = i+1 ; j4 = j+1
               cycle 
            case(2)
               i2 = i ; j2 = j+1
               i3 = i-1 ; j3 = j
               i4 = i ; j4 = j+1
            case(3)
               i2 = i-1 ; j2 = j
               i3 = i ; j3 = j-1
               i4 = i ; j4 = j    
            case(4)
               ! i2 = i ; j2 = j-1
               ! i3 = i+1 ; j3 = j
               ! i4 = i+1 ; j4 = j 
               cycle
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
            Gama_r = Gama_r + 1.d0/sound(i4,j4)/den(i4,j4)
            right_hand = right_hand + 1.d0/sound(i4,j4)/den(i4,j4)*pre(i4,j4) - &
                        (velx(i4,j4)*(cos(phi)) + vely(i4,j4)*(sin(phi))) / sin(delphi/2.d0)
         enddo       
         Pstar1(i,j) = right_hand / Gama_r     

      Case ('given_velocity')      
      
         Vstarx(i,j) = tvel(1)  
         Vstary(i,j) = tvel(2)
         
         right_hand  = 0.d0
         Gama_r      = 0.d0         
         do side = 1 , 4      !> 4 means quadrangle
            select case(side)
            case(1)
               ! i2 = i+1 ; j2 = j
               ! i3 = i ; j3 = j + 1
               ! i4 = i+1 ; j4 = j+1
               cycle 
            case(2)
               i2 = i ; j2 = j+1
               i3 = i-1 ; j3 = j
               i4 = i ; j4 = j+1
            case(3)
               i2 = i-1 ; j2 = j
               i3 = i ; j3 = j-1
               i4 = i ; j4 = j    
            case(4)
               ! i2 = i ; j2 = j-1
               ! i3 = i+1 ; j3 = j
               ! i4 = i+1 ; j4 = j 
               cycle
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
            Gama_r = Gama_r + 1.d0/sound(i4,j4)/den(i4,j4)
            right_hand = right_hand + 1.d0/sound(i4,j4)/den(i4,j4)*pre(i4,j4) - &
                        (velx(i4,j4)*(cos(phi)) + vely(i4,j4)*(sin(phi))) / sin(delphi/2.d0)
         enddo       
         Pstar1(i,j) = right_hand / Gama_r       
      
      Case ('given_pressure')
         
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
               ! i2 = i+1 ; j2 = j
               ! i3 = i ; j3 = j + 1
               ! i4 = i+1 ; j4 = j+1
               cycle 
            case(2)
               i2 = i ; j2 = j+1
               i3 = i-1 ; j3 = j
               i4 = i ; j4 = j+1
            case(3)
               i2 = i-1 ; j2 = j
               i3 = i ; j3 = j-1
               i4 = i ; j4 = j    
            case(4)
               ! i2 = i ; j2 = j-1
               ! i3 = i+1 ; j3 = j
               ! i4 = i+1 ; j4 = j 
               cycle
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
                  
            SMx = SMx + Vjr11*velx(i4,j4) + Vjr12*vely(i4,j4) - &
                  (pre(i4,j4)-tpre)*cos(phi)*sin(delphi/2.d0)
            SMy = SMy + Vjr12*velx(i4,j4) + Vjr22*vely(i4,j4) - &
                  (pre(i4,j4)-tpre)*sin(phi)*sin(delphi/2.d0) 
            
         enddo
         Vstarx(i,j)=(SMx*temptB-SMy*temptC)/(temptA*temptB-temptC*temptC) 
         Vstary(i,j)=(SMy*temptA-SMx*temptC)/(temptA*temptB-temptC*temptC)      
         Pstar1(i,j) = tpre !right_hand / Gama_r
         
         vx = -vertx(i,j) / sqrt(vertx(i,j)**2 + verty(i,j)**2)
         vy = -verty(i,j) / sqrt(vertx(i,j)**2 + verty(i,j)**2)


         theta_m = acos(vx/sqrt(vx*vx+vy*vy))
         right_hand = 0.d0
         Gama_r     = 0.d0  
    
         do side = 1 , 4      !> 4 means quadrangle
            select case(side)
            case(1)
               ! i2 = i+1 ; j2 = j
               ! i3 = i ; j3 = j + 1
               ! i4 = i+1 ; j4 = j+1
               cycle 
            case(2)
               i2 = i ; j2 = j+1
               i3 = i-1 ; j3 = j
               i4 = i ; j4 = j+1
            case(3)
               i2 = i-1 ; j2 = j
               i3 = i ; j3 = j-1
               i4 = i ; j4 = j    
            case(4)
               ! i2 = i ; j2 = j-1
               ! i3 = i+1 ; j3 = j
               ! i4 = i+1 ; j4 = j 
               cycle
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
            Gama_r = Gama_r + (vx / sqrt(vx*vx+vy*vy)*cos(phi) + vy / sqrt(vx*vx+vy*vy)*sin(phi))/cos(delphi/2.d0)
            right_hand = right_hand + (pre(i4,j4)-tpre-Pstar1(i,j))/sound(i4,j4)/den(i4,j4) - &
                        (velx(i4,j4)*(cos(phi)) + vely(i4,j4)*(sin(phi)))/sin(delphi/2.d0)
         enddo            
         u_bar = right_hand / Gama_r       
         Vstarx(i,j) = abs(u_bar) * vx / sqrt(vx*vx+vy*vy)
         Vstary(i,j) = abs(u_bar) * vy / sqrt(vx*vx+vy*vy)        

      Case ('slide_line')
      
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
               ! i2 = i+1 ; j2 = j
               ! i3 = i ; j3 = j + 1
               ! i4 = i+1 ; j4 = j+1
               cycle 
            case(2)
               i2 = i ; j2 = j+1
               i3 = i-1 ; j3 = j
               i4 = i ; j4 = j+1
            case(3)
               i2 = i-1 ; j2 = j
               i3 = i ; j3 = j-1
               i4 = i ; j4 = j    
            case(4)
               ! i2 = i ; j2 = j-1
               ! i3 = i+1 ; j3 = j
               ! i4 = i+1 ; j4 = j 
               cycle
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

         A(1,1) = temptA
         A(1,2) = temptC
         A(2,1) = temptC
         A(2,2) = temptB
         B(1)   = SMx
         B(2)   = SMy
         call tangent_otimes_tangent(tangent,tot)
         call tot_dot_A_dot_tot(tot,A,tot_A_tot)
         barA = tot_A_tot + (temptA+temptB)*(Ione - tot)            
         call tot_dot_B(tot,B,barB)
         vx = (barB(1)*barA(2,2)-barB(2)*barA(1,2))/(barA(1,1)*barA(2,2)-barA(1,2)*barA(1,2))
         vy = (barB(2)*barA(1,1)-barB(1)*barA(1,2))/(barA(1,1)*barA(2,2)-barA(1,2)*barA(1,2))
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
               ! i2 = i+1 ; j2 = j
               ! i3 = i ; j3 = j + 1
               ! i4 = i+1 ; j4 = j+1
               cycle 
            case(2)
               i2 = i ; j2 = j+1
               i3 = i-1 ; j3 = j
               i4 = i ; j4 = j+1
            case(3)
               i2 = i-1 ; j2 = j
               i3 = i ; j3 = j-1
               i4 = i ; j4 = j    
            case(4)
               ! i2 = i ; j2 = j-1
               ! i3 = i+1 ; j3 = j
               ! i4 = i+1 ; j4 = j 
               cycle
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
      Case default
         
         print*, 'node_type_char(i,j)=',node_type_char(i,j)
         stop 'node_type_char(i,j) is error, pls call xihua f8a'
   
      End Select
      
   Enddo

   !> y_bgn
   j = 0
   Do i = 1 , nx - 1
      tangent = node_type_direction(i,j,1:2)
      tvel    = node_type_velocity(i,j,1:2)
      tpre    = node_type_pressure(i,j)    
      Select case(node_type_char(i,j))
      
      Case ('fix')
         Vstarx(i,j) = 0.d0  
         Vstary(i,j) = 0.d0 
         
         right_hand  = 0.d0
         Gama_r      = 0.d0         
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
               ! i2 = i-1 ; j2 = j
               ! i3 = i ; j3 = j-1
               ! i4 = i ; j4 = j     
               cycle
            case(4)
               ! i2 = i ; j2 = j-1
               ! i3 = i+1 ; j3 = j
               ! i4 = i+1 ; j4 = j  
               cycle
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
            Gama_r = Gama_r + 1.d0/sound(i4,j4)/den(i4,j4)
            right_hand = right_hand + 1.d0/sound(i4,j4)/den(i4,j4)*pre(i4,j4) - &
                        (velx(i4,j4)*(cos(phi)) + vely(i4,j4)*(sin(phi))) / sin(delphi/2.d0)
         enddo       
         Pstar1(i,j) = right_hand / Gama_r     

      Case ('given_velocity')      
         Vstarx(i,j) = tvel(1)  
         Vstary(i,j) = tvel(2)
         
         right_hand  = 0.d0
         Gama_r      = 0.d0         
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
               ! i2 = i-1 ; j2 = j
               ! i3 = i ; j3 = j-1
               ! i4 = i ; j4 = j     
               cycle
            case(4)
               ! i2 = i ; j2 = j-1
               ! i3 = i+1 ; j3 = j
               ! i4 = i+1 ; j4 = j  
               cycle
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
            Gama_r = Gama_r + 1.d0/sound(i4,j4)/den(i4,j4)
            right_hand = right_hand + 1.d0/sound(i4,j4)/den(i4,j4)*pre(i4,j4) - &
                        (velx(i4,j4)*(cos(phi)) + vely(i4,j4)*(sin(phi))) / sin(delphi/2.d0)
         enddo       
         Pstar1(i,j) = right_hand / Gama_r       
      
      Case ('given_pressure')      
         stop 'Not Reay, 1154fdai'
      Case ('slide_line')
      
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
               ! i2 = i-1 ; j2 = j
               ! i3 = i ; j3 = j-1
               ! i4 = i ; j4 = j     
               cycle
            case(4)
               ! i2 = i ; j2 = j-1
               ! i3 = i+1 ; j3 = j
               ! i4 = i+1 ; j4 = j  
               cycle
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

         A(1,1) = temptA
         A(1,2) = temptC
         A(2,1) = temptC
         A(2,2) = temptB
         B(1)   = SMx
         B(2)   = SMy
         call tangent_otimes_tangent(tangent,tot)
         call tot_dot_A_dot_tot(tot,A,tot_A_tot)
         barA = tot_A_tot + (temptA+temptB)*(Ione - tot)            
         call tot_dot_B(tot,B,barB)
         vx = (barB(1)*barA(2,2)-barB(2)*barA(1,2))/(barA(1,1)*barA(2,2)-barA(1,2)*barA(1,2))
         vy = (barB(2)*barA(1,1)-barB(1)*barA(1,2))/(barA(1,1)*barA(2,2)-barA(1,2)*barA(1,2))
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
               ! i2 = i-1 ; j2 = j
               ! i3 = i ; j3 = j-1
               ! i4 = i ; j4 = j     
               cycle
            case(4)
               ! i2 = i ; j2 = j-1
               ! i3 = i+1 ; j3 = j
               ! i4 = i+1 ; j4 = j  
               cycle
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
      Case default
         
         print*, 'node_type_char(i,j)=',node_type_char(i,j)
         stop 'node_type_char(i,j) is error, pls call xihua f8a'
   
      End Select
      
   Enddo
   
   !> y_end
   j = ny
   Do i = 1 , nx - 1
      tangent = node_type_direction(i,j,1:2)
      tvel    = node_type_velocity(i,j,1:2)
      tpre    = node_type_pressure(i,j)    
      Select case(node_type_char(i,j))
      
      Case ('fix')
         Vstarx(i,j) = 0.d0  
         Vstary(i,j) = 0.d0 
         
         right_hand  = 0.d0
         Gama_r      = 0.d0         
         do side = 1 , 4      !> 4 means quadrangle
            select case(side)
            case(1)
               ! i2 = i+1 ; j2 = j
               ! i3 = i ; j3 = j + 1
               ! i4 = i+1 ; j4 = j+1 
               cycle
            case(2)
               ! i2 = i ; j2 = j+1
               ! i3 = i-1 ; j3 = j
               ! i4 = i ; j4 = j+1
               cycle 
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
            Gama_r = Gama_r + 1.d0/sound(i4,j4)/den(i4,j4)
            right_hand = right_hand + 1.d0/sound(i4,j4)/den(i4,j4)*pre(i4,j4) - &
                        (velx(i4,j4)*(cos(phi)) + vely(i4,j4)*(sin(phi))) / sin(delphi/2.d0)
         enddo       
         Pstar1(i,j) = right_hand / Gama_r     

      Case ('given_velocity')      
         Vstarx(i,j) = tvel(1)  
         Vstary(i,j) = tvel(2)
         
         right_hand  = 0.d0
         Gama_r      = 0.d0         
         do side = 1 , 4      !> 4 means quadrangle
            select case(side)
            case(1)
               ! i2 = i+1 ; j2 = j
               ! i3 = i ; j3 = j + 1
               ! i4 = i+1 ; j4 = j+1 
               cycle
            case(2)
               ! i2 = i ; j2 = j+1
               ! i3 = i-1 ; j3 = j
               ! i4 = i ; j4 = j+1
               cycle 
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
            Gama_r = Gama_r + 1.d0/sound(i4,j4)/den(i4,j4)
            right_hand = right_hand + 1.d0/sound(i4,j4)/den(i4,j4)*pre(i4,j4) - &
                        (velx(i4,j4)*(cos(phi)) + vely(i4,j4)*(sin(phi))) / sin(delphi/2.d0)
         enddo       
         Pstar1(i,j) = right_hand / Gama_r       
      
      Case ('given_pressure')      
         stop 'Not Ready xihua1428'
      Case ('slide_line')
      
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
               ! i2 = i+1 ; j2 = j
               ! i3 = i ; j3 = j + 1
               ! i4 = i+1 ; j4 = j+1 
               cycle
            case(2)
               ! i2 = i ; j2 = j+1
               ! i3 = i-1 ; j3 = j
               ! i4 = i ; j4 = j+1
               cycle 
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

         A(1,1) = temptA
         A(1,2) = temptC
         A(2,1) = temptC
         A(2,2) = temptB
         B(1)   = SMx
         B(2)   = SMy
         call tangent_otimes_tangent(tangent,tot)
         call tot_dot_A_dot_tot(tot,A,tot_A_tot)
         barA = tot_A_tot + (temptA+temptB)*(Ione - tot)            
         call tot_dot_B(tot,B,barB)
         vx = (barB(1)*barA(2,2)-barB(2)*barA(1,2))/(barA(1,1)*barA(2,2)-barA(1,2)*barA(1,2))
         vy = (barB(2)*barA(1,1)-barB(1)*barA(1,2))/(barA(1,1)*barA(2,2)-barA(1,2)*barA(1,2))
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
               ! i2 = i+1 ; j2 = j
               ! i3 = i ; j3 = j + 1
               ! i4 = i+1 ; j4 = j+1 
               cycle
            case(2)
               ! i2 = i ; j2 = j+1
               ! i3 = i-1 ; j3 = j
               ! i4 = i ; j4 = j+1
               cycle 
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
      Case default
         
         print*, 'node_type_char(i,j)=',node_type_char(i,j)
         stop 'node_type_char(i,j) is error, pls call xihua f8a'
   
      End Select
      
   Enddo

   !> x_bgn & y_bgn
   i = 0
   Do j = 0 , 0
      tangent = node_type_direction(i,j,1:2)
      tvel    = node_type_velocity(i,j,1:2)
      tpre    = node_type_pressure(i,j)   

      if(Gou_Xing == 'Quadrangle' .or. Origin_YesNo == 'No') then
      
         Select case(node_type_char(i,j))
         
         Case ('fix')
            Vstarx(i,j) = 0.d0  
            Vstary(i,j) = 0.d0 
            
            right_hand  = 0.d0
            Gama_r      = 0.d0         
            do side = 1 , 4      !> 4 means quadrangle
               select case(side)
               case(1)
                  i2 = i+1 ; j2 = j
                  i3 = i ; j3 = j + 1
                  i4 = i+1 ; j4 = j+1
               case(2)
                  ! i2 = i ; j2 = j+1
                  ! i3 = i-1 ; j3 = j
                  ! i4 = i ; j4 = j+1
                  cycle
               case(3)
                  ! i2 = i-1 ; j2 = j
                  ! i3 = i ; j3 = j-1
                  ! i4 = i ; j4 = j     
                  cycle
               case(4)
                  ! i2 = i ; j2 = j-1
                  ! i3 = i+1 ; j3 = j
                  ! i4 = i+1 ; j4 = j 
                  cycle               
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
               Gama_r = Gama_r + 1.d0/sound(i4,j4)/den(i4,j4)
               right_hand = right_hand + 1.d0/sound(i4,j4)/den(i4,j4)*pre(i4,j4) - &
                           (velx(i4,j4)*(cos(phi)) + vely(i4,j4)*(sin(phi))) / sin(delphi/2.d0)
            enddo       
            Pstar1(i,j) = right_hand / Gama_r     

         Case ('given_velocity')      
            Vstarx(i,j) = tvel(1)  
            Vstary(i,j) = tvel(2)
            
            right_hand  = 0.d0
            Gama_r      = 0.d0         
            do side = 1 , 4      !> 4 means quadrangle
               select case(side)
               case(1)
                  i2 = i+1 ; j2 = j
                  i3 = i ; j3 = j + 1
                  i4 = i+1 ; j4 = j+1
               case(2)
                  ! i2 = i ; j2 = j+1
                  ! i3 = i-1 ; j3 = j
                  ! i4 = i ; j4 = j+1
                  cycle
               case(3)
                  ! i2 = i-1 ; j2 = j
                  ! i3 = i ; j3 = j-1
                  ! i4 = i ; j4 = j     
                  cycle
               case(4)
                  ! i2 = i ; j2 = j-1
                  ! i3 = i+1 ; j3 = j
                  ! i4 = i+1 ; j4 = j 
                  cycle               
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
               Gama_r = Gama_r + 1.d0/sound(i4,j4)/den(i4,j4)
               right_hand = right_hand + 1.d0/sound(i4,j4)/den(i4,j4)*pre(i4,j4) - &
                           (velx(i4,j4)*(cos(phi)) + vely(i4,j4)*(sin(phi))) / sin(delphi/2.d0)
            enddo       
            Pstar1(i,j) = right_hand / Gama_r       
         
         Case ('given_pressure')
            stop 'not ready xihua1707'         
        
         Case default
            
            print*, 'node_type_char(i,j)=',node_type_char(i,j)
            stop 'node_type_char(i,j) is error, pls call xihua f8a'
      
         End Select
         
      elseif(Gou_Xing == 'Sector' .and. Origin_YesNo == 'Yes') then
      
         Select case(node_type_char(i,j))
         
         Case ('fix')
            Vstarx(i,j) = 0.d0  
            Vstary(i,j) = 0.d0 
            
            right_hand  = 0.d0
            Gama_r      = 0.d0         
            do side = 1 , 1      !> 4 means quadrangle
               i2 = i+1 ; j2 = j
               i3 = i+1 ; j3 = j + 1
               i4 = i+1 ; j4 = j+1     
         
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
               Gama_r = Gama_r + 1.d0/sound(i4,j4)/den(i4,j4)
               right_hand = right_hand + 1.d0/sound(i4,j4)/den(i4,j4)*pre(i4,j4) - &
                           (velx(i4,j4)*(cos(phi)) + vely(i4,j4)*(sin(phi))) !/ sin(delphi/2.d0)
            enddo       
            Pstar1(i,j) = right_hand / Gama_r     

         Case ('given_velocity')      
            Vstarx(i,j) = tvel(1)  
            Vstary(i,j) = tvel(2)
            
            right_hand  = 0.d0
            Gama_r      = 0.d0         
            do side = 1 , 1      !> 4 means quadrangle
               i2 = i+1 ; j2 = j
               i3 = i+1 ; j3 = j + 1
               i4 = i+1 ; j4 = j+1 
         
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
               Gama_r = Gama_r + 1.d0/sound(i4,j4)/den(i4,j4)
               right_hand = right_hand + 1.d0/sound(i4,j4)/den(i4,j4)*pre(i4,j4) - &
                           (velx(i4,j4)*(cos(phi)) + vely(i4,j4)*(sin(phi))) / sin(delphi/2.d0)
            enddo       
            Pstar1(i,j) = right_hand / Gama_r       
         
         Case ('given_pressure')
            stop 'not ready xihua1707'         
         
         Case default
            
            print*, 'node_type_char(i,j)=',node_type_char(i,j)
            stop 'node_type_char(i,j) is error, pls call xihua f8a'
      
         End Select
               
      else 
         stop 'xihuafdfda9a'
         
      endif 
   Enddo
   
   !> x_bgn & y_end
   i = 0
   Do j = ny, ny 
      tangent = node_type_direction(i,j,1:2)
      tvel    = node_type_velocity(i,j,1:2)
      tpre    = node_type_pressure(i,j)   

      if(Gou_Xing == 'Quadrangle' .or. Origin_YesNo == 'No') then
      
         Select case(node_type_char(i,j))
         
         Case ('fix')
            Vstarx(i,j) = 0.d0  
            Vstary(i,j) = 0.d0 
            
            right_hand  = 0.d0
            Gama_r      = 0.d0         
            do side = 1 , 4      !> 4 means quadrangle
               select case(side)
               case(1)
                  ! i2 = i+1 ; j2 = j
                  ! i3 = i ; j3 = j + 1
                  ! i4 = i+1 ; j4 = j+1
                  cycle 
               case(2)
                  ! i2 = i ; j2 = j+1
                  ! i3 = i-1 ; j3 = j
                  ! i4 = i ; j4 = j+1
                  cycle
               case(3)
                  ! i2 = i-1 ; j2 = j
                  ! i3 = i ; j3 = j-1
                  ! i4 = i ; j4 = j     
                  cycle
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
               Gama_r = Gama_r + 1.d0/sound(i4,j4)/den(i4,j4)
               right_hand = right_hand + 1.d0/sound(i4,j4)/den(i4,j4)*pre(i4,j4) - &
                           (velx(i4,j4)*(cos(phi)) + vely(i4,j4)*(sin(phi))) / sin(delphi/2.d0)
            enddo       
            Pstar1(i,j) = right_hand / Gama_r     

         Case ('given_velocity')      
            Vstarx(i,j) = tvel(1)  
            Vstary(i,j) = tvel(2)
            
            right_hand  = 0.d0
            Gama_r      = 0.d0         
            do side = 1 , 4      !> 4 means quadrangle
               select case(side)
               case(1)
                  ! i2 = i+1 ; j2 = j
                  ! i3 = i ; j3 = j + 1
                  ! i4 = i+1 ; j4 = j+1
                  cycle 
               case(2)
                  ! i2 = i ; j2 = j+1
                  ! i3 = i-1 ; j3 = j
                  ! i4 = i ; j4 = j+1
                  cycle
               case(3)
                  ! i2 = i-1 ; j2 = j
                  ! i3 = i ; j3 = j-1
                  ! i4 = i ; j4 = j     
                  cycle
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
               Gama_r = Gama_r + 1.d0/sound(i4,j4)/den(i4,j4)
               right_hand = right_hand + 1.d0/sound(i4,j4)/den(i4,j4)*pre(i4,j4) - &
                           (velx(i4,j4)*(cos(phi)) + vely(i4,j4)*(sin(phi))) / sin(delphi/2.d0)
            enddo       
            Pstar1(i,j) = right_hand / Gama_r       
         
         Case ('given_pressure')
            stop 'not ready xihua1707'         
        
         Case default
            
            print*, 'node_type_char(i,j)=',node_type_char(i,j)
            stop 'node_type_char(i,j) is error, pls call xihua f8a'
      
         End Select
   
      elseif(Gou_Xing == 'Sector' .and. Origin_YesNo == 'Yes')  then
      
         Select case(node_type_char(i,j))
         
         Case ('fix')
            Vstarx(i,j) = 0.d0  
            Vstary(i,j) = 0.d0 
            
            right_hand  = 0.d0
            Gama_r      = 0.d0         
            do side = 1 , 4      !> 4 means quadrangle
               select case(side)
               case(1)
                  ! i2 = i+1 ; j2 = j
                  ! i3 = i ; j3 = j + 1
                  ! i4 = i+1 ; j4 = j+1
                  cycle 
               case(2)
                  ! i2 = i ; j2 = j+1
                  ! i3 = i-1 ; j3 = j
                  ! i4 = i ; j4 = j+1
                  cycle
               case(3)
                  ! i2 = i-1 ; j2 = j
                  ! i3 = i ; j3 = j-1
                  ! i4 = i ; j4 = j     
                  cycle
               case(4)
                  i2 = i+1 ; j2 = j-1
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
               Gama_r = Gama_r + 1.d0/sound(i4,j4)/den(i4,j4)
               right_hand = right_hand + 1.d0/sound(i4,j4)/den(i4,j4)*pre(i4,j4) - &
                           (velx(i4,j4)*(cos(phi)) + vely(i4,j4)*(sin(phi))) !/ sin(delphi/2.d0)
            enddo       
            Pstar1(i,j) = right_hand / Gama_r     

         Case ('given_velocity')      
            Vstarx(i,j) = tvel(1)  
            Vstary(i,j) = tvel(2)
            
            right_hand  = 0.d0
            Gama_r      = 0.d0         
            do side = 1 , 4      !> 4 means quadrangle
               select case(side)
               case(1)
                  ! i2 = i+1 ; j2 = j
                  ! i3 = i ; j3 = j + 1
                  ! i4 = i+1 ; j4 = j+1
                  cycle 
               case(2)
                  ! i2 = i ; j2 = j+1
                  ! i3 = i-1 ; j3 = j
                  ! i4 = i ; j4 = j+1
                  cycle
               case(3)
                  ! i2 = i-1 ; j2 = j
                  ! i3 = i ; j3 = j-1
                  ! i4 = i ; j4 = j     
                  cycle
               case(4)
                  i2 = i+1 ; j2 = j-1
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
               Gama_r = Gama_r + 1.d0/sound(i4,j4)/den(i4,j4)
               right_hand = right_hand + 1.d0/sound(i4,j4)/den(i4,j4)*pre(i4,j4) - &
                           (velx(i4,j4)*(cos(phi)) + vely(i4,j4)*(sin(phi))) / sin(delphi/2.d0)
            enddo       
            Pstar1(i,j) = right_hand / Gama_r       
         
         Case ('given_pressure')
            stop 'not ready xihua1707'         
         
         Case default
            
            print*, 'node_type_char(i,j)=',node_type_char(i,j)
            stop 'node_type_char(i,j) is error, pls call xihua f8a'
      
         End Select
               
      else 
         stop 'xihuafdfda9a'
         
      endif 
   Enddo
   
   !> x_end & y_bgn
   i = nx
   Do j = 0 , 0
      tangent = node_type_direction(i,j,1:2)
      tvel    = node_type_velocity(i,j,1:2)
      tpre    = node_type_pressure(i,j)    
      Select case(node_type_char(i,j))
      
      Case ('fix')
         Vstarx(i,j) = 0.d0  
         Vstary(i,j) = 0.d0 
         
         right_hand  = 0.d0
         Gama_r      = 0.d0         
         do side = 1 , 4      !> 4 means quadrangle
            select case(side)
            case(1)
               ! i2 = i+1 ; j2 = j
               ! i3 = i ; j3 = j + 1
               ! i4 = i+1 ; j4 = j+1
               cycle 
            case(2)
               i2 = i ; j2 = j+1
               i3 = i-1 ; j3 = j
               i4 = i ; j4 = j+1
            case(3)
               ! i2 = i-1 ; j2 = j
               ! i3 = i ; j3 = j-1
               ! i4 = i ; j4 = j    
               cycle
            case(4)
               ! i2 = i ; j2 = j-1
               ! i3 = i+1 ; j3 = j
               ! i4 = i+1 ; j4 = j 
               cycle
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
            Gama_r = Gama_r + 1.d0/sound(i4,j4)/den(i4,j4)
            right_hand = right_hand + 1.d0/sound(i4,j4)/den(i4,j4)*pre(i4,j4) - &
                        (velx(i4,j4)*(cos(phi)) + vely(i4,j4)*(sin(phi))) / sin(delphi/2.d0)
         enddo       
         Pstar1(i,j) = right_hand / Gama_r     

      Case ('given_velocity')      
         Vstarx(i,j) = tvel(1)  
         Vstary(i,j) = tvel(2)
         
         right_hand  = 0.d0
         Gama_r      = 0.d0         
         do side = 1 , 4      !> 4 means quadrangle
            select case(side)
            case(1)
               ! i2 = i+1 ; j2 = j
               ! i3 = i ; j3 = j + 1
               ! i4 = i+1 ; j4 = j+1
               cycle 
            case(2)
               i2 = i ; j2 = j+1
               i3 = i-1 ; j3 = j
               i4 = i ; j4 = j+1
            case(3)
               ! i2 = i-1 ; j2 = j
               ! i3 = i ; j3 = j-1
               ! i4 = i ; j4 = j    
               cycle
            case(4)
               ! i2 = i ; j2 = j-1
               ! i3 = i+1 ; j3 = j
               ! i4 = i+1 ; j4 = j 
               cycle
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
            Gama_r = Gama_r + 1.d0/sound(i4,j4)/den(i4,j4)
            right_hand = right_hand + 1.d0/sound(i4,j4)/den(i4,j4)*pre(i4,j4) - &
                        (velx(i4,j4)*(cos(phi)) + vely(i4,j4)*(sin(phi))) / sin(delphi/2.d0)
         enddo       
         Pstar1(i,j) = right_hand / Gama_r       
      
      Case ('given_pressure')
            stop 'not ready xihua1707'         
      
      case('mix_slide_line_given_pressure') 
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
               ! i2 = i+1 ; j2 = j
               ! i3 = i ; j3 = j + 1
               ! i4 = i+1 ; j4 = j+1
               cycle 
            case(2)
               i2 = i ; j2 = j+1
               i3 = i-1 ; j3 = j
               i4 = i ; j4 = j+1
            case(3)
               ! i2 = i-1 ; j2 = j
               ! i3 = i ; j3 = j-1
               ! i4 = i ; j4 = j    
               cycle
            case(4)
               ! i2 = i ; j2 = j-1
               ! i3 = i+1 ; j3 = j
               ! i4 = i+1 ; j4 = j 
               cycle
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
                  
            SMx = SMx + Vjr11*velx(i4,j4) + Vjr12*vely(i4,j4) - &
                  (pre(i4,j4)-tpre)*cos(phi)*sin(delphi/2.d0)
            SMy = SMy + Vjr12*velx(i4,j4) + Vjr22*vely(i4,j4) - &
                  (pre(i4,j4)-tpre)*sin(phi)*sin(delphi/2.d0)  
            
            phi = phi2 + delphi/2.d0
            Gama_r = Gama_r + 1.d0/sound(i4,j4)/den(i4,j4)
            right_hand = right_hand + 1.d0/sound(i4,j4)/den(i4,j4)*(pre(i4,j4)-tpre) - &
                        (velx(i4,j4)*(cos(phi)) + vely(i4,j4)*(sin(phi))) / sin(delphi/2.d0)           

         enddo
         Vstarx(i,j)=(SMx*temptB-SMy*temptC)/(temptA*temptB-temptC*temptC) 
         Vstary(i,j)=(SMy*temptA-SMx*temptC)/(temptA*temptB-temptC*temptC)      
         Pstar1(i,j) = tpre     

         A(1,1) = temptA
         A(1,2) = temptC
         A(2,1) = temptC
         A(2,2) = temptB
         B(1)   = SMx
         B(2)   = SMy
         call tangent_otimes_tangent(tangent,tot)
         call tot_dot_A_dot_tot(tot,A,tot_A_tot)
         barA = tot_A_tot + (temptA+temptB)*(Ione - tot)            
         call tot_dot_B(tot,B,barB)
         vx = (barB(1)*barA(2,2)-barB(2)*barA(1,2))/(barA(1,1)*barA(2,2)-barA(1,2)*barA(1,2))
         vy = (barB(2)*barA(1,1)-barB(1)*barA(1,2))/(barA(1,1)*barA(2,2)-barA(1,2)*barA(1,2))
         if( sqrt(vx*vx+vy*vy) < 1.d-10) then
            Vstarx(i,j) = 0.d0
            Vstary(i,j) = 0.d0 
            cycle
         endif 
         Vstarx(i,j) = vx
         Vstary(i,j) = vy
         
         theta_m = acos(vx/sqrt(vx*vx+vy*vy))       
         
         right_hand = 0.d0
         Gama_r     = 0.d0  
    
         do side = 1 , 4      !> 4 means quadrangle
            select case(side)
            case(1)
               ! i2 = i+1 ; j2 = j
               ! i3 = i ; j3 = j + 1
               ! i4 = i+1 ; j4 = j+1
               cycle 
            case(2)
               i2 = i ; j2 = j+1
               i3 = i-1 ; j3 = j
               i4 = i ; j4 = j+1
            case(3)
               ! i2 = i-1 ; j2 = j
               ! i3 = i ; j3 = j-1
               ! i4 = i ; j4 = j    
               cycle
            case(4)
               ! i2 = i ; j2 = j-1
               ! i3 = i+1 ; j3 = j
               ! i4 = i+1 ; j4 = j 
               cycle
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
            Gama_r = Gama_r + (vx / sqrt(vx*vx+vy*vy)*cos(phi) + vy / sqrt(vx*vx+vy*vy)*sin(phi))/cos(delphi/2.d0)
            right_hand = right_hand + (pre(i4,j4)-tpre-Pstar1(i,j))/sound(i4,j4)/den(i4,j4) - &
                        (velx(i4,j4)*(cos(phi)) + vely(i4,j4)*(sin(phi)))/sin(delphi/2.d0)
         enddo            
         u_bar = right_hand / Gama_r       
         Vstarx(i,j) = abs(u_bar) * vx / sqrt(vx*vx+vy*vy)
         Vstary(i,j) = abs(u_bar) * vy / sqrt(vx*vx+vy*vy)          

      Case default
         
         print*, 'node_type_char(i,j)=',node_type_char(i,j)
         stop 'node_type_char(i,j) is error, pls call xihua f8a'
   
      End Select
      
   Enddo

   !> x_end & y_end
   i = nx
   Do j = ny , ny 
      tangent = node_type_direction(i,j,1:2)
      tvel    = node_type_velocity(i,j,1:2)
      tpre    = node_type_pressure(i,j)    
      Select case(node_type_char(i,j))
      
      Case ('fix')
         Vstarx(i,j) = 0.d0  
         Vstary(i,j) = 0.d0 
         
         right_hand  = 0.d0
         Gama_r      = 0.d0         
         do side = 1 , 4      !> 4 means quadrangle
            select case(side)
            case(1)
               ! i2 = i+1 ; j2 = j
               ! i3 = i ; j3 = j + 1
               ! i4 = i+1 ; j4 = j+1
               cycle 
            case(2)
               ! i2 = i ; j2 = j+1
               ! i3 = i-1 ; j3 = j
               ! i4 = i ; j4 = j+1
               cycle 
            case(3)
               i2 = i-1 ; j2 = j
               i3 = i ; j3 = j-1
               i4 = i ; j4 = j    
            case(4)
               ! i2 = i ; j2 = j-1
               ! i3 = i+1 ; j3 = j
               ! i4 = i+1 ; j4 = j 
               cycle
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
            Gama_r = Gama_r + 1.d0/sound(i4,j4)/den(i4,j4)
            right_hand = right_hand + 1.d0/sound(i4,j4)/den(i4,j4)*pre(i4,j4) - &
                        (velx(i4,j4)*(cos(phi)) + vely(i4,j4)*(sin(phi))) / sin(delphi/2.d0)
         enddo       
         Pstar1(i,j) = right_hand / Gama_r     

      Case ('given_velocity')      
         Vstarx(i,j) = tvel(1)  
         Vstary(i,j) = tvel(2)
         
         right_hand  = 0.d0
         Gama_r      = 0.d0         
         do side = 1 , 4      !> 4 means quadrangle
            select case(side)
            case(1)
               ! i2 = i+1 ; j2 = j
               ! i3 = i ; j3 = j + 1
               ! i4 = i+1 ; j4 = j+1
               cycle 
            case(2)
               ! i2 = i ; j2 = j+1
               ! i3 = i-1 ; j3 = j
               ! i4 = i ; j4 = j+1
               cycle 
            case(3)
               i2 = i-1 ; j2 = j
               i3 = i ; j3 = j-1
               i4 = i ; j4 = j    
            case(4)
               ! i2 = i ; j2 = j-1
               ! i3 = i+1 ; j3 = j
               ! i4 = i+1 ; j4 = j 
               cycle
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
            Gama_r = Gama_r + 1.d0/sound(i4,j4)/den(i4,j4)
            right_hand = right_hand + 1.d0/sound(i4,j4)/den(i4,j4)*pre(i4,j4) - &
                        (velx(i4,j4)*(cos(phi)) + vely(i4,j4)*(sin(phi))) / sin(delphi/2.d0)
         enddo       
         Pstar1(i,j) = right_hand / Gama_r       
      
      Case ('given_pressure')
            stop 'not ready xihuafewe'               
      
      Case ('mix_slide_line_given_pressure')
      
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
               ! i2 = i+1 ; j2 = j
               ! i3 = i ; j3 = j + 1
               ! i4 = i+1 ; j4 = j+1
               cycle 
            case(2)
               ! i2 = i ; j2 = j+1
               ! i3 = i-1 ; j3 = j
               ! i4 = i ; j4 = j+1
               cycle 
            case(3)
               i2 = i-1 ; j2 = j
               i3 = i ; j3 = j-1
               i4 = i ; j4 = j    
            case(4)
               ! i2 = i ; j2 = j-1
               ! i3 = i+1 ; j3 = j
               ! i4 = i+1 ; j4 = j 
               cycle
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
                  
            SMx = SMx + Vjr11*velx(i4,j4) + Vjr12*vely(i4,j4) - &
                  (pre(i4,j4)-tpre)*cos(phi)*sin(delphi/2.d0)
            SMy = SMy + Vjr12*velx(i4,j4) + Vjr22*vely(i4,j4) - &
                  (pre(i4,j4)-tpre)*sin(phi)*sin(delphi/2.d0)  
            
            phi = phi2 + delphi/2.d0
            Gama_r = Gama_r + 1.d0/sound(i4,j4)/den(i4,j4)
            right_hand = right_hand + 1.d0/sound(i4,j4)/den(i4,j4)*(pre(i4,j4)-tpre) - &
                        (velx(i4,j4)*(cos(phi)) + vely(i4,j4)*(sin(phi))) / sin(delphi/2.d0)           

         enddo
         Vstarx(i,j)=(SMx*temptB-SMy*temptC)/(temptA*temptB-temptC*temptC) 
         Vstary(i,j)=(SMy*temptA-SMx*temptC)/(temptA*temptB-temptC*temptC)      
         Pstar1(i,j) = tpre   

         A(1,1) = temptA
         A(1,2) = temptC
         A(2,1) = temptC
         A(2,2) = temptB
         B(1)   = SMx
         B(2)   = SMy
         call tangent_otimes_tangent(tangent,tot)
         call tot_dot_A_dot_tot(tot,A,tot_A_tot)
         barA = tot_A_tot + (temptA+temptB)*(Ione - tot)            
         call tot_dot_B(tot,B,barB)
         vx = (barB(1)*barA(2,2)-barB(2)*barA(1,2))/(barA(1,1)*barA(2,2)-barA(1,2)*barA(1,2))
         vy = (barB(2)*barA(1,1)-barB(1)*barA(1,2))/(barA(1,1)*barA(2,2)-barA(1,2)*barA(1,2))
         if( sqrt(vx*vx+vy*vy) < 1.d-10) then
            Vstarx(i,j) = 0.d0
            Vstary(i,j) = 0.d0 
            cycle
         endif 
         Vstarx(i,j) = vx
         Vstary(i,j) = vy
         
         theta_m = acos(vx/sqrt(vx*vx+vy*vy)) 
                  
         right_hand = 0.d0
         Gama_r     = 0.d0  
    
         do side = 1 , 4      !> 4 means quadrangle
            select case(side)
            case(1)
               ! i2 = i+1 ; j2 = j
               ! i3 = i ; j3 = j + 1
               ! i4 = i+1 ; j4 = j+1
               cycle 
            case(2)
               ! i2 = i ; j2 = j+1
               ! i3 = i-1 ; j3 = j
               ! i4 = i ; j4 = j+1
               cycle 
            case(3)
               i2 = i-1 ; j2 = j
               i3 = i ; j3 = j-1
               i4 = i ; j4 = j    
            case(4)
               ! i2 = i ; j2 = j-1
               ! i3 = i+1 ; j3 = j
               ! i4 = i+1 ; j4 = j 
               cycle
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
            Gama_r = Gama_r + (vx / sqrt(vx*vx+vy*vy)*cos(phi) + vy / sqrt(vx*vx+vy*vy)*sin(phi))/cos(delphi/2.d0)
            right_hand = right_hand + (pre(i4,j4)-tpre-Pstar1(i,j))/sound(i4,j4)/den(i4,j4) - &
                        (velx(i4,j4)*(cos(phi)) + vely(i4,j4)*(sin(phi)))/sin(delphi/2.d0)
         enddo            
         u_bar = right_hand / Gama_r       
         Vstarx(i,j) = abs(u_bar) * vx / sqrt(vx*vx+vy*vy)
         Vstary(i,j) = abs(u_bar) * vy / sqrt(vx*vx+vy*vy)  
         
      Case default
         
         print*, 'node_type_char(i,j)=',node_type_char(i,j)
         stop 'node_type_char(i,j) is error, pls call xihua f8a'
   
      End Select
      
   Enddo

End Subroutine Node_Solver_Cat_UeAS_BndyNode_add_pian

