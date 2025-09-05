
!> 20250701 
!> xihua
!> 
subroutine cal_geo_2d(nx,ny,vertx,verty,den,area,volm,subcell_area_ca,&
   smass,amass,parea,pvolm,pmass,smass_tri)
   implicit none
   integer  :: nx,ny 
   real*8,intent(in)    :: vertx(0:nx,0:ny),verty(0:nx,0:ny),den(nx,ny)
   real*8,intent(out)   :: area(nx,ny),volm(nx,ny),subcell_area_ca(nx,ny,4)
   real*8,intent(out)   :: smass(nx,ny,4),amass(nx,ny),smass_tri(nx,ny,8)
   real*8,intent(out)   :: parea(0:nx,0:ny),pvolm(0:nx,0:ny),pmass(0:nx,0:ny)

   real*8  :: subcell_area_ca_tri(nx,ny,8)
   integer  :: i,j 
   real*8   :: xr1,yr1,xr2,yr2,ca_area,cy_volm

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

      call cal_subcell_area_ca(vertx(i-1,j-1),verty(i-1,j-1), &
            vertx(i,j-1),verty(i,j-1),vertx(i,j),verty(i,j),  & 
            vertx(i-1,j),verty(i-1,j),subcell_area_ca(i,j,:), &
            subcell_area_ca_tri(i,j,:))  
      smass(i,j,:) = den(i,j) * subcell_area_ca(i,j,:)
      smass_tri(i,j,:) = den(i,j) * subcell_area_ca_tri(i,j,:)
      amass(i,j)= sum(smass(i,j,:))
	enddo
	enddo

   do i = 1 , nx - 1
   do j = 1 , ny - 1 
      pmass(i,j) = smass(i,j,3) + smass(i+1,j,4) &
                 + smass(i+1,j+1,1) + smass(i,j+1,2)
      pvolm(i,j) = subcell_area_ca(i,j,3)     + subcell_area_ca(i+1,j,4) &
                 + subcell_area_ca(i+1,j+1,1) + subcell_area_ca(i,j+1,2) 
      parea(i,j) = subcell_area_ca(i,j,3)     + subcell_area_ca(i+1,j,4) &
                 + subcell_area_ca(i+1,j+1,1) + subcell_area_ca(i,j+1,2) 
   enddo 
   enddo 

   do i = 0 , 0
   do j = 1 , ny - 1 
      pmass(i,j) = smass(i+1,j,4) + smass(i+1,j+1,1) 
      parea(i,j) = subcell_area_ca(i+1,j,4) + subcell_area_ca(i+1,j+1,1)
      pvolm(i,j) = parea(i,j)  
   enddo 
   enddo 

   do i = nx , nx 
   do j = 1 , ny - 1 
      pmass(i,j) = smass(i,j,3) + smass(i,j+1,2)
      parea(i,j) = subcell_area_ca(i,j,3)     + subcell_area_ca(i,j+1,2) 
      pvolm(i,j) = parea(i,j) 
   enddo 
   enddo 

   do i = 1 , nx - 1
   do j = 0 , 0 
      pmass(i,j) = smass(i+1,j+1,1) + smass(i,j+1,2)
      parea(i,j) = subcell_area_ca(i+1,j+1,1) + subcell_area_ca(i,j+1,2) 
      pvolm(i,j) = parea(i,j) 
   enddo 
   enddo 

   do i = 1 , nx - 1
   do j = ny , ny 
      pmass(i,j) = smass(i,j,3) + smass(i+1,j,4)  
      parea(i,j) = subcell_area_ca(i,j,3)     + subcell_area_ca(i+1,j,4)                
      pvolm(i,j) = parea(i,j) 
   enddo 
   enddo 

   do i = 0 , 0
   do j = 0 , 0  
      pmass(i,j) = smass(i+1,j+1,1) 
      parea(i,j) = subcell_area_ca(i+1,j+1,1)  
      pvolm(i,j) = parea(i,j) 
   enddo 
   enddo 

   do i = 0,0
   do j = ny,ny 
      pmass(i,j) = smass(i+1,j,4) 
      parea(i,j) = subcell_area_ca(i+1,j,4) 
      pvolm(i,j) = parea(i,j) 
   enddo 
   enddo 

   do i = nx,nx
   do j = 0,0
      pmass(i,j) = smass(i,j+1,2)
      parea(i,j) = subcell_area_ca(i,j+1,2) 
      pvolm(i,j) = parea(i,j) 
   enddo 
   enddo 

   do i = nx,nx 
   do j = ny,ny  
      pmass(i,j) = smass(i,j,3)
      parea(i,j) = subcell_area_ca(i,j,3)
      pvolm(i,j) = parea(i,j) 
   enddo 
   enddo 

end subroutine cal_geo_2d


subroutine triangle_geo_f1(x1,y1,x2,y2,xc,yc,A12,A2c,Ac1)
   implicit none
   real*8,intent(in)    :: x1,y1,x2,y2,xc,yc
   real*8,intent(out)   :: A12(2),A2c(2),Ac1(2)

   real*8   :: alength12,anormalx12,anormaly12
   real*8   :: alength2c,anormalx2c,anormaly2c
   real*8   :: alengthc1,anormalxc1,anormalyc1,area

   call cal_length_normal(x1,x2,y1,y2,alength12,anormalx12,anormaly12)
   call cal_length_normal(x2,xc,y2,yc,alength2c,anormalx2c,anormaly2c)
   call cal_length_normal(xc,x1,yc,y1,alengthc1,anormalxc1,anormalyc1)
   
   A12(1) = alengthc1 * anormalxc1 + alength2c * anormalx2c
   A12(2) = alengthc1 * anormalyc1 + alength2c * anormaly2c

   A2c(1) = alengthc1 * anormalxc1 + alength12 * anormalx12
   A2c(2) = alengthc1 * anormalyc1 + alength12 * anormaly12

   Ac1(1) = alength2c * anormalx2c + alength12 * anormalx12
   Ac1(2) = alength2c * anormaly2c + alength12 * anormaly12

   A12 =  A12 / alength12
   A2c =  A2c / alength12
   Ac1 =  Ac1 / alength12

end subroutine triangle_geo_f1

subroutine triangle_geo_f2(xx1,yy1,xx2,yy2,xxc,yyc,uu1,vv1,uu2,vv2,BBB)
   implicit none
   real*8,intent(in)    :: xx1,yy1,xx2,yy2,xxc,yyc,uu1,vv1,uu2,vv2
   real*8,intent(out)   :: BBB

   real*8   :: xx5,yy5,alength15,anormalx15,anormaly15
   real*8   :: alength52,anormalx52,anormaly52
   real*8   :: alengthc5,anormalxc5,anormalyc5
   real*8   :: alength5c,anormalx5c,anormaly5c

   xx5 = ( xx1 + xx2 ) / 2.d0
   yy5 = ( yy1 + yy2 ) / 2.d0

   call cal_length_normal(xx1,xx5,yy1,yy5,alength15,anormalx15,anormaly15)
   call cal_length_normal(xx5,xx2,yy5,yy2,alength52,anormalx52,anormaly52)
   call cal_length_normal(xxc,xx5,yyc,yy5,alengthc5,anormalxc5,anormalyc5)
   call cal_length_normal(xx5,xxc,yy5,yyc,alength5c,anormalx5c,anormaly5c)

   BBB = alength15*(anormalx15*uu1 + anormaly15*vv1) + &
         alength5c*(anormalx5c*uu1 + anormaly5c*vv1) + &
         alength52*(anormalx52*uu2 + anormaly52*vv2) + &
         alengthc5*(anormalxc5*uu2 + anormalyc5*vv2)

end subroutine triangle_geo_f2



subroutine cal_length_normal(x1,x2,y1,y2,alength,anormalx,anormaly)
	implicit none
	real*8 :: x1,x2,y1,y2
	real*8 :: alength,anormalx,anormaly

	call Length(x1,y1,x2,y2,alength)
   if ( alength < 1.d-14 ) then 
      anormalx = 0.d0 
      anormaly = 0.d0 
   else 
	   call Normal(x1,y1,x2,y2,anormalx,anormaly)
	   anormalx = -anormalx
	   anormaly = -anormaly		
   endif 

end subroutine cal_length_normal

!> 20170923
!> xihua 
!> compute the length of segment

subroutine Length(x1,y1,x2,y2,alength)
   implicit none
   real*8,intent(in) :: x1,y1,x2,y2
   real*8,intent(out) :: alength

   alength = dsqrt((x1-x2)**2+(y1-y2)**2)    
   
end subroutine Length

!> 20170923
!> xihua 
!> compute the unit outward normal of segment

subroutine Normal(x1,y1,x2,y2,anormalx,anormaly)
   implicit none
   real*8,intent(in) :: x1,y1,x2,y2
   real*8,intent(out) :: anormalx,anormaly
   real*8 :: tempt

   tempt=(x1-x2)**2+(y1-y2)**2
   if( tempt < 1.d-14 ) then 
      anormalx = 0.d0 
      anormaly = 0.d0 
   else 
      tempt = sqrt(tempt)
      anormalx=(y1-y2)/tempt
      anormaly=(x2-x1)/tempt
   endif 
end subroutine Normal


subroutine cal_eos(eos,gamar)
   implicit none 
   integer,intent(in)  :: eos 
   real*8,intent(out)  :: gamar

   select case (eos)
   case (715)
      gamar = 7.d0 / 5.d0 
   case (513)
      gamar = 5.d0 / 3.d0 
   case (815)
      gamar = 8.d0 / 5.d0 
   case (312)
      gamar = 3.d0 / 2.d0
   case default
      stop 'eos,errorela'
   end select 
end subroutine cal_eos
            
Real*8 Function ca_area(xa,ya,xb,yb)
   implicit none
   real*8,intent(in)  :: xa,ya
   real*8,intent(in)  :: xb,yb
   ca_area = (xa-xb)*(ya + yb)/2.d0
End Function ca_area

! Real*8 Function ca_volm(xa,ya,xb,yb)
!    implicit none
!    real*8,intent(in)  :: xa,ya
!    real*8,intent(in)  :: xb,yb
!    ca_volm = (xa*yb - xb*ya)/2.d0
! End Function ca_volm

subroutine ca_area_quad(x1,y1,x2,y2,x3,y3,x4,y4,area)
   implicit none
   real*8,intent(in)  :: x1,y1,x2,y2,x3,y3,x4,y4
   real*8,intent(out) :: area 
   real*8  :: ca_area
   area = 0.d0 
   area = area + ca_area(x1,y1,x2,y2) & 
               + ca_area(x2,y2,x3,y3) &
               + ca_area(x3,y3,x4,y4) &
               + ca_area(x4,y4,x1,y1) 
end subroutine ca_area_quad


subroutine ca_area_triangle(x1,y1,x2,y2,x3,y3,area)
   implicit none
   real*8,intent(in)  :: x1,y1,x2,y2,x3,y3
   real*8,intent(out) :: area 
   real*8  :: ca_area
   area = 0.d0 
   area = area + ca_area(x1,y1,x2,y2) & 
               + ca_area(x2,y2,x3,y3) &
               + ca_area(x3,y3,x1,y1) 
end subroutine ca_area_triangle

subroutine cal_subcell_area_ca(x1,y1,x2,y2,x3,y3,x4,y4,subcell_area,subcell_area_tri)
   implicit none 
   real*8,intent(in)    :: x1,y1,x2,y2,x3,y3,x4,y4
   real*8,intent(out)   :: subcell_area(4),subcell_area_tri(8)
   real*8   :: xc,yc,x12,y12,x23,y23,x34,y34,x41,y41
   subcell_area = 0.d0 
   xc = (x1+x2+x3+x4)/4.d0 
   yc = (y1+y2+y3+y4)/4.d0 
   x12 = (x1+x2)/2.d0 
   y12 = (y1+y2)/2.d0 
   x23 = (x3+x2)/2.d0 
   y23 = (y3+y2)/2.d0 
   x34 = (x3+x4)/2.d0 
   y34 = (y3+y4)/2.d0 
   x41 = (x1+x4)/2.d0 
   y41 = (y1+y4)/2.d0
   
   call cal_quad_cell_area_ca(x1,y1,x12,y12,xc,yc,x41,y41,subcell_area(1)) 
   call cal_quad_cell_area_ca(x2,y2,x23,y23,xc,yc,x12,y12,subcell_area(2)) 
   call cal_quad_cell_area_ca(x3,y3,x34,y34,xc,yc,x23,y23,subcell_area(3)) 
   call cal_quad_cell_area_ca(x4,y4,x41,y41,xc,yc,x34,y34,subcell_area(4)) 

   call cal_tri_cell_area_ca(x1,y1,x12,y12,xc,yc,subcell_area_tri(1)) 
   call cal_tri_cell_area_ca(x12,y12,x2,y2,xc,yc,subcell_area_tri(2))   
   call cal_tri_cell_area_ca(x2,y2,x23,y23,xc,yc,subcell_area_tri(3)) 
   call cal_tri_cell_area_ca(x23,y23,x3,y3,xc,yc,subcell_area_tri(4)) 
   call cal_tri_cell_area_ca(x3,y3,x34,y34,xc,yc,subcell_area_tri(5)) 
   call cal_tri_cell_area_ca(x34,y34,x4,y4,xc,yc,subcell_area_tri(6))
   call cal_tri_cell_area_ca(x4,y4,x41,y41,xc,yc,subcell_area_tri(7)) 
   call cal_tri_cell_area_ca(x41,y41,x1,y1,xc,yc,subcell_area_tri(8)) 

end subroutine cal_subcell_area_ca

subroutine cal_tri_cell_area_ca(x1,y1,x2,y2,x3,y3,area)
   implicit none 
   real*8,intent(in)    :: x1,y1,x2,y2,x3,y3
   real*8,intent(out)   :: area 
   real*8  :: ca_area
   area = 0.d0 
   area = area + ca_area(x1,y1,x2,y2) & 
               + ca_area(x2,y2,x3,y3) &
               + ca_area(x3,y3,x1,y1) 
end subroutine cal_tri_cell_area_ca


subroutine cal_quad_cell_area_ca(x1,y1,x2,y2,x3,y3,x4,y4,area)
   implicit none 
   real*8,intent(in)    :: x1,y1,x2,y2,x3,y3,x4,y4
   real*8,intent(out)   :: area 
   real*8  :: ca_area
   area = 0.d0 
   area = area + ca_area(x1,y1,x2,y2) & 
               + ca_area(x2,y2,x3,y3) &
               + ca_area(x3,y3,x4,y4) &
               + ca_area(x4,y4,x1,y1) 
end subroutine cal_quad_cell_area_ca

Real*8 Function cy_volm(za,ra,zb,rb)
   implicit none
   real*8,intent(in)  :: za,ra
   real*8,intent(in)  :: zb,rb
   cy_volm = (za-zb)*(ra*ra + ra*rb + rb*rb)/6.d0
End Function cy_volm

subroutine cal_subcell_volm_cy(x1,y1,x2,y2,x3,y3,x4,y4,subcell_volm)
   implicit none 
   real*8,intent(in)    :: x1,y1,x2,y2,x3,y3,x4,y4
   real*8,intent(out)   :: subcell_volm(4)
   real*8   :: xc,yc,x12,y12,x23,y23,x34,y34,x41,y41
   subcell_volm = 0.d0 
   xc = (x1+x2+x3+x4)/4.d0 
   yc = (y1+y2+y3+y4)/4.d0 
   x12 = (x1+x2)/2.d0 
   y12 = (y1+y2)/2.d0 
   x23 = (x3+x2)/2.d0 
   y23 = (y3+y2)/2.d0 
   x34 = (x3+x4)/2.d0 
   y34 = (y3+y4)/2.d0 
   x41 = (x1+x4)/2.d0 
   y41 = (y1+y4)/2.d0
   call cal_quad_cell_volm_cy(x1,y1,x12,y12,xc,yc,x41,y41,subcell_volm(1)) 
   call cal_quad_cell_volm_cy(x2,y2,x23,y23,xc,yc,x12,y12,subcell_volm(2)) 
   call cal_quad_cell_volm_cy(x3,y3,x34,y34,xc,yc,x23,y23,subcell_volm(3)) 
   call cal_quad_cell_volm_cy(x4,y4,x41,y41,xc,yc,x34,y34,subcell_volm(4)) 
end subroutine cal_subcell_volm_cy

subroutine cal_quad_cell_volm_cy(x1,y1,x2,y2,x3,y3,x4,y4,volm)
   implicit none 
   real*8,intent(in)    :: x1,y1,x2,y2,x3,y3,x4,y4
   real*8,intent(out)   :: volm 
   real*8  :: cy_volm
   volm = 0.d0 
   volm = volm + cy_volm(x1,y1,x2,y2) & 
               + cy_volm(x2,y2,x3,y3) &
               + cy_volm(x3,y3,x4,y4) &
               + cy_volm(x4,y4,x1,y1) 
end subroutine cal_quad_cell_volm_cy


subroutine cal_geo_quad_vel(x1,y1,x2,y2,x3,y3,x4,y4,u1,v1,u2,v2,u4,v4, &
                        u3,v3,sum_vel_ln)
   implicit none 
	real*8,intent(in)   :: x1,y1,x2,y2,x3,y3,x4,y4
	real*8,intent(in)   :: u1,v1,u2,v2,u3,v3,u4,v4
	real*8,intent(out)  :: sum_vel_ln

   real*8 :: ln12(2),ln23(2),ln34(2),ln41(2)
   real*8 :: alength,anormalx,anormaly

   call cal_length_normal(x1,x2,y1,y2,alength,anormalx,anormaly)
   ln12(1) = alength * anormalx
   ln12(2) = alength * anormaly
   
   call cal_length_normal(x2,x3,y2,y3,alength,anormalx,anormaly)
   ln23(1) = alength * anormalx
   ln23(2) = alength * anormaly

   call cal_length_normal(x3,x4,y3,y4,alength,anormalx,anormaly)
   ln34(1) = alength * anormalx
   ln34(2) = alength * anormaly

   call cal_length_normal(x4,x1,y4,y1,alength,anormalx,anormaly)
   ln41(1) = alength * anormalx
   ln41(2) = alength * anormaly

   sum_vel_ln = ln12(1)*(u1+u2) + ln12(2)*(v1+v2) + &
                ln23(1)*(u2+u3) + ln23(2)*(v2+v3) + &
                ln34(1)*(u3+u4) + ln34(2)*(v3+v4) + &
                ln41(1)*(u4+u1) + ln41(2)*(v4+v1) 
   sum_vel_ln = sum_vel_ln / 2.d0 

end subroutine cal_geo_quad_vel

subroutine cal_geo_quad(x1,y1,x2,y2,x3,y3,x4,y4, &
            u1,v1,u2,v2,u4,v4,hbar,area,tmpB,vec214,vec234)
   implicit none 
	real*8,intent(in)   :: x1,y1,x2,y2,x3,y3,x4,y4
   real*8,intent(in)   :: u1,v1,u2,v2,u4,v4
	real*8,intent(out)  :: hbar,area,vec214(2),vec234(2),tmpB

   real*8 :: ln12(2),ln23(2),ln34(2),ln41(2)
   real*8 :: alength,anormalx,anormaly

   hbar = ( sqrt( (x3-x1)**2 + (y3-y1)**2 ) + &
            sqrt( (x2-x4)**2 + (y2-y4)**2 )  ) / 2.d0
   area = 0.d0 

   call cal_length_normal(x1,x2,y1,y2,alength,anormalx,anormaly)
   ln12(1) = alength * anormalx
   ln12(2) = alength * anormaly
   area = area + 0.5d0*(x1*y2-x2*y1)

   call cal_length_normal(x2,x3,y2,y3,alength,anormalx,anormaly)
   ln23(1) = alength * anormalx
   ln23(2) = alength * anormaly
   area = area + 0.5d0*(x2*y3-x3*y2)

   call cal_length_normal(x3,x4,y3,y4,alength,anormalx,anormaly)
   ln34(1) = alength * anormalx
   ln34(2) = alength * anormaly
   area = area + 0.5d0*(x3*y4-x4*y3)

   call cal_length_normal(x4,x1,y4,y1,alength,anormalx,anormaly)
   ln41(1) = alength * anormalx
   ln41(2) = alength * anormaly
   area = area + 0.5d0*(x4*y1-x1*y4)

   tmpB = u1*(ln12(1)+ln41(1)) + v1*(ln12(2)+ln41(2)) + &
          u2*(ln12(1)+ln23(1)) + v2*(ln12(2)+ln23(2)) + &
          u4*(ln34(1)+ln41(1)) + v4*(ln41(2)+ln41(2)) 

   vec214 = ln12 + ln41 
   vec234 = ln23 + ln34

end subroutine cal_geo_quad

subroutine cal_mass_center(x1,y1,x2,y2,x3,y3,x4,y4,xc,yc)
	implicit none 
	real*8,intent(in)   :: x1,y1,x2,y2,x3,y3,x4,y4
   real*8,intent(out)  :: xc,yc 

   real*8   :: lc,area 

   call cal_quad_cell_area_ca(x1,y1,x2,y2,x3,y3,x4,y4,area)

   xc = 0.d0 
   call cal_int_l_x_dxy(x1,y1,x2,y2,lc) ; xc = xc + lc
   call cal_int_l_x_dxy(x2,y2,x3,y3,lc) ; xc = xc + lc
   call cal_int_l_x_dxy(x3,y3,x4,y4,lc) ; xc = xc + lc
   call cal_int_l_x_dxy(x4,y4,x1,y1,lc) ; xc = xc + lc
   xc = xc / area

   yc = 0.d0 
   call cal_int_l_y_dxy(x1,y1,x2,y2,lc) ; yc = yc + lc
   call cal_int_l_y_dxy(x2,y2,x3,y3,lc) ; yc = yc + lc
   call cal_int_l_y_dxy(x3,y3,x4,y4,lc) ; yc = yc + lc
   call cal_int_l_y_dxy(x4,y4,x1,y1,lc) ; yc = yc + lc   
   yc = yc / area 

end subroutine cal_mass_center

subroutine cal_int_l_x_dxy(x1,y1,x2,y2,lx)
   implicit none 
   real*8,intent(in)    :: x1,y1,x2,y2 
   real*8,intent(out)   :: lx

   lx = (y2-y1)*(x1*x1+x1*x2+x2*x2)/6.d0

end subroutine cal_int_l_x_dxy

subroutine cal_int_l_y_dxy(x1,y1,x2,y2,ly)
   implicit none 
   real*8,intent(in)    :: x1,y1,x2,y2 
   real*8,intent(out)   :: ly

   ly = (x1-x2)*(y1*y1+y1*y2+y2*y2)/6.d0

end subroutine cal_int_l_y_dxy