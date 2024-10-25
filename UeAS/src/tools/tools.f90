
Subroutine Cal_Geometry(Ini_CDT,MeshCY2D)
   use COMMON_TYPE_Condition
   use COMMON_TYPE_Mesh
   implicit none
   
   Type(COMMON_TYPE_CDT),intent(in) :: Ini_CDT
   Type(Mesh_Cyl_2D),intent(inout) :: MeshCY2D
   
   integer :: ix, iy   
   real*8 :: xr1,xr2,yr1,yr2,cyl_area,cyl_volm,cat_area,cat_volm

   Select Case(Ini_CDT%Cat_Cyl)

   Case('Cat')

      do ix = 1 , Ini_CDT%nx
      do iy = 1 , Ini_CDT%ny
         MeshCY2D%area(ix,iy) = 0.d0
         MeshCY2D%volm(ix,iy) = 0.d0

         xr1=MeshCY2D%vertex_x(ix-1, iy-1)
         yr1=MeshCY2D%vertex_y(ix-1, iy-1)
         xr2=MeshCY2D%vertex_x(ix  , iy-1)
         yr2=MeshCY2D%vertex_y(ix  , iy-1)
         MeshCY2D%area(ix,iy) = MeshCY2D%area(ix,iy) + cat_area(xr1,yr1,xr2,yr2)
         MeshCY2D%volm(ix,iy) = MeshCY2D%volm(ix,iy) + cat_volm(xr1,yr1,xr2,yr2)   
         
         xr1=MeshCY2D%vertex_x(ix  , iy-1)
         yr1=MeshCY2D%vertex_y(ix  , iy-1)
         xr2=MeshCY2D%vertex_x(ix  , iy)
         yr2=MeshCY2D%vertex_y(ix  , iy)
         MeshCY2D%area(ix,iy) = MeshCY2D%area(ix,iy) + cat_area(xr1,yr1,xr2,yr2)
         MeshCY2D%volm(ix,iy) = MeshCY2D%volm(ix,iy) + cat_volm(xr1,yr1,xr2,yr2)  
         
         xr1=MeshCY2D%vertex_x(ix  , iy)
         yr1=MeshCY2D%vertex_y(ix  , iy)
         xr2=MeshCY2D%vertex_x(ix-1, iy)
         yr2=MeshCY2D%vertex_y(ix-1, iy)
         MeshCY2D%area(ix,iy) = MeshCY2D%area(ix,iy) + cat_area(xr1,yr1,xr2,yr2)
         MeshCY2D%volm(ix,iy) = MeshCY2D%volm(ix,iy) + cat_volm(xr1,yr1,xr2,yr2)     

         xr1=MeshCY2D%vertex_x(ix-1, iy)
         yr1=MeshCY2D%vertex_y(ix-1, iy)
         xr2=MeshCY2D%vertex_x(ix-1, iy-1)
         yr2=MeshCY2D%vertex_y(ix-1, iy-1)
         MeshCY2D%area(ix,iy) = MeshCY2D%area(ix,iy) + cat_area(xr1,yr1,xr2,yr2)
         MeshCY2D%volm(ix,iy) = MeshCY2D%volm(ix,iy) + cat_volm(xr1,yr1,xr2,yr2)

      enddo
      enddo   
      
   Case default
      
      stop 'Cat or Cyl is Error, pls Call Xihua.'
     
   EndSelect          
      

End Subroutine Cal_Geometry

   
Real*8 Function cat_area(xa,ra,xb,rb)
   implicit none
   real*8,intent(in)  :: xa,ra,xb,rb

   cat_area = (xa*rb - xb*ra)/2.d0

End Function cat_area

Real*8 Function cat_volm(xa,ra,xb,rb)
   implicit none
   real*8,intent(in)  :: xa,ra,xb,rb
   
   cat_volm = (xa*rb - xb*ra)/2.d0

End Function cat_volm

Real*8 Function cyl_area(xa,ra,xb,rb)
   implicit none
   real*8,intent(in)  :: xa,ra,xb,rb

   cyl_area = (xa-xb)*(ra + rb)/2.d0

End Function cyl_area

Real*8 Function cyl_volm(xa,ra,xb,rb)
   implicit none
   real*8,intent(in)  :: xa,ra,xb,rb
   
   cyl_volm = (xa-xb)*(ra*ra + ra*rb + rb*rb)/6.d0

End Function cyl_volm


subroutine Length(x1,y1,x2,y2,alength)
   implicit none 
   real*8,intent(in)  :: x1,y1,x2,y2
   real*8,intent(out) :: alength
   real*8 :: r1,r2

   alength =  dsqrt((x1-x2)**2+(y1-y2)**2) 

end subroutine Length
      
      
subroutine Normal(x1,y1,x2,y2,anormalx,anormaly)
   implicit none 
   real*8,intent(in)  :: x1,y1,x2,y2
   real*8,intent(out) :: anormalx,anormaly

   anormalx=(y1-y2)/dsqrt((x1-x2)**2+(y1-y2)**2)
   anormaly=(x2-x1)/dsqrt((x1-x2)**2+(y1-y2)**2)

end subroutine Normal


subroutine three(a11,a12,a13,a21,a22,a23,a31,a32,a33,det)
   implicit none
   real*8,intent(in)  :: a11,a12,a13,a21,a22,a23,a31,a32,a33
   real*8,intent(out) :: det

   det = a11*(a22*a33-a23*a32) - &
         a12*(a21*a33-a23*a31) + &
         a13*(a21*a32-a22*a31)

end subroutine three
   
	
subroutine Cal_theta(x,y,theta)
	implicit none
	real*8,intent(in) :: x,y
	real*8,intent(out) :: theta
	real*8 :: eps,PI,a,b,r
	
	eps = 1.d-30
	PI = 3.1415926535897932384626d0
	if(dabs(x)<eps) then
		if(dabs(y)<eps) then
			stop 'X,Y is (0,0)'
		elseif(y>eps) then
			theta = PI/2.d0
		elseif(y<-eps) then
			theta = 3.d0*PI/2.d0
		endif
		return
	endif
	
	if(dabs(y)<eps) then
		if(x>eps) then
			theta = 0.d0
		elseif(x<-eps) then
			theta = PI
		endif
		return
	endif
	
	r = dsqrt(x*x + y*y)
	a = dacos(x/r)
	b = dasin(y/r)
	if(a>eps .and. b>eps) then      
		theta = a  
	elseif(a>eps .and. b<-eps) then  
		theta = 2.d0*PI - a
	elseif(a<-eps .and. b>eps) then  
		theta = PI - a
	elseif(a<-eps .and. b>eps) then  
		theta = PI + a
	endif
	
End subroutine Cal_theta
	
subroutine Cal_njr(phi1,phi2,njr1,njr2)
	implicit none
	real*8,intent(in) :: phi1,phi2
	real*8,intent(out) :: njr1,njr2
	
	njr1 = -dsin(phi2) + dsin(phi1)
	njr2 =  dcos(phi2) - dcos(phi1)
End Subroutine Cal_njr


subroutine jifen_theta2co2stheta(x1,x2,y)
   implicit none
   real*8,intent(in)  :: x1,x2
   real*8,intent(out) :: y
   real*8 :: y1,y2
   
   y2= 1.d0/6.d0*x2**3 + 1.d0/4.d0*x2**2*dsin(2.d0*x2) + &
       1.d0/4.d0*x2*dcos(2.d0*x2) - 1.d0/8.d0*dsin(2.d0*x2)

   y1= 1.d0/6.d0*x1**3 + 1.d0/4.d0*x1**2*dsin(2.d0*x1) + &
       1.d0/4.d0*x1*dcos(2.d0*x1) - 1.d0/8.d0*dsin(2.d0*x1)
       
   y = y2 - y1
   
end subroutine jifen_theta2co2stheta


subroutine jifen_theta2costhetasintheta(x1,x2,y)
   implicit none
   real*8,intent(in)  :: x1,x2
   real*8,intent(out) :: y
   real*8 :: y1,y2
   
   y2= -1.d0/4.d0*x2**2*cos(2.d0*x2) + 1.d0/4.d0*x2*sin(2.d0*x2) + &
        1.d0/8.d0*cos(2.d0*x2)

   y1= -1.d0/4.d0*x1**2*cos(2.d0*x1) + 1.d0/4.d0*x1*sin(2.d0*x1) + &
        1.d0/8.d0*cos(2.d0*x1)
       
   y = y2 - y1
   
end subroutine jifen_theta2costhetasintheta

subroutine jifen_theta2si2ntheta(x1,x2,y)
   implicit none
   real*8,intent(in)  :: x1,x2
   real*8,intent(out) :: y
   real*8 :: y1,y2
   
   y2= 1.d0/6.d0*x2**3 - 1.d0/4.d0*x2**2*sin(2.d0*x2) - &
       1.d0/4.d0*x2*cos(2.d0*x2) + 1.d0/8.d0*sin(2.d0*x2)

   y1= 1.d0/6.d0*x1**3 - 1.d0/4.d0*x1**2*sin(2.d0*x1) - &
       1.d0/4.d0*x1*cos(2.d0*x1) + 1.d0/8.d0*sin(2.d0*x1)
       
   y = y2 - y1
   
end subroutine jifen_theta2si2ntheta

subroutine jifen_thetaco2stheta(x1,x2,y)
   implicit none
   real*8,intent(in)  :: x1,x2
   real*8,intent(out) :: y
   real*8 :: y1,y2
   
   y2= (2.d0*x2*sin(2.d0*x2) + cos(2.d0*x2) + 2.d0*x2**2 )/8.d0

   y1= (2.d0*x1*sin(2.d0*x1) + cos(2.d0*x1) + 2.d0*x1**2 )/8.d0

   y = y2 - y1
   
end subroutine jifen_thetaco2stheta

subroutine jifen_thetasi2ntheta(x1,x2,y)
   implicit none
   real*8,intent(in)  :: x1,x2
   real*8,intent(out) :: y
   real*8 :: y1,y2
   
   y2= x2**2/4.d0 - x2*sin(2.d0*x2)/4.d0 - cos(2.d0*x2)/8.d0

   y1= x1**2/4.d0 - x1*sin(2.d0*x1)/4.d0 - cos(2.d0*x1)/8.d0

   y = y2 - y1
   
end subroutine jifen_thetasi2ntheta

subroutine jifen_thetacosthetasintheta(x1,x2,y)
   implicit none
   real*8,intent(in)  :: x1,x2
   real*8,intent(out) :: y
   real*8 :: y1,y2
   
   y2= -x2*cos(2.d0*x2)/4.d0 + sin(2.d0*x2)/8.d0

   y1= -x1*cos(2.d0*x1)/4.d0 + sin(2.d0*x1)/8.d0

   y = y2 - y1
   
end subroutine jifen_thetacosthetasintheta

subroutine jifen_co2stheta(x1,x2,y)
   implicit none
   real*8,intent(in)  :: x1,x2
   real*8,intent(out) :: y
   real*8 :: y1,y2
   
   y2= x2/2.d0 + sin(2.d0*x2)/4.d0

   y1= x1/2.d0 + sin(2.d0*x1)/4.d0

   y = y2 - y1
   
end subroutine jifen_co2stheta

subroutine jifen_si2ntheta(x1,x2,y)
   implicit none
   real*8,intent(in)  :: x1,x2
   real*8,intent(out) :: y
   real*8 :: y1,y2
   
   y2= x2/2.d0 - sin(2.d0*x2)/4.d0

   y1= x1/2.d0 - sin(2.d0*x1)/4.d0

   y = y2 - y1
   
end subroutine jifen_si2ntheta

subroutine jifen_costhetasintheta(x1,x2,y)
   implicit none
   real*8,intent(in)  :: x1,x2
   real*8,intent(out) :: y
   real*8 :: y1,y2
   
   y2= -cos(2.d0*x2)/4.d0

   y1= -cos(2.d0*x1)/4.d0

   y = y2 - y1
   
end subroutine jifen_costhetasintheta

subroutine jifen_thetacostheta(x1,x2,y)
   implicit none
   real*8,intent(in)  :: x1,x2
   real*8,intent(out) :: y
   real*8 :: y1,y2
   
   y2= x2*sin(x2) + cos(x2)

   y1= x1*sin(x1) + cos(x1)

   y = y2 - y1
   
end subroutine jifen_thetacostheta


subroutine jifen_thetasintheta(x1,x2,y)
   implicit none
   real*8,intent(in)  :: x1,x2
   real*8,intent(out) :: y
   real*8 :: y1,y2
   
   y2= -x2*cos(x2) + sin(x2)

   y1= -x1*cos(x1) + sin(x1)

   y = y2 - y1
   
end subroutine jifen_thetasintheta

   
Subroutine a_otimes_b(a,b,c)
   implicit none 
   real*8,intent(in)  :: a(2),b(2)
   real*8,intent(out) :: c(2,2)
   c(1,1) = a(1)*a(1)
   c(1,2) = a(1)*b(2)
   c(2,1) = a(2)*b(1)
   c(2,2) = b(2)*b(2)
End Subroutine a_otimes_b

Subroutine matrix_times_matrix(a,b,c)
   implicit none
   real*8,intent(in)  :: a(2,2),b(2,2)
   real*8,intent(out) :: c(2,2)
   c(1,1) = a(1,1)*b(1,1) + a(1,2)*b(2,1)
   c(1,2) = a(1,1)*b(2,1) + a(1,2)*b(2,2)
   c(2,1) = a(2,1)*b(1,1) + a(2,2)*b(2,1)
   c(2,2) = a(2,1)*b(2,1) + a(2,2)*b(2,2)   
End Subroutine matrix_times_matrix

Subroutine matrix_times_vec(a,b,c)
   implicit none
   real*8,intent(in)  :: a(2,2),b(2)
   real*8,intent(out) :: c(2)
   c(1) = a(1,1)*b(1) + a(1,2)*b(2)
   c(2) = a(2,1)*b(1) + a(2,2)*b(2)   
End Subroutine matrix_times_vec   
     
subroutine tangent_otimes_tangent(tangent,tot)
   implicit none 
   real*8,intent(in)  :: tangent(2)
   real*8,intent(out) :: tot(2,2)
   
   call a_otimes_b(tangent,tangent,tot)
   
end subroutine tangent_otimes_tangent

subroutine tot_dot_A_dot_tot(tot,A,tot_A_tot)     
  implicit none 
  real*8,intent(in)   :: tot(2,2),A(2,2)
  real*8,intent(out)  :: tot_A_tot(2,2)
  
  real*8  :: tmp(2,2)
  
  call matrix_times_matrix(tot,A,tmp)
  call matrix_times_matrix(tmp,tot,tot_A_tot)
  
end subroutine tot_dot_A_dot_tot


subroutine tot_dot_B(tot,B,barB)
   implicit none 
   real*8,intent(in)   :: tot(2,2),B(2)
   real*8,intent(out)  :: barB(2)
   
   call matrix_times_vec(tot,B,barB)

end subroutine tot_dot_B
   

