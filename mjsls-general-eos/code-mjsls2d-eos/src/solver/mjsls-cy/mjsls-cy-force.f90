

!> 20230112
!> xihua
!> compute the force on the nodal 
   
subroutine cal_node_force_cy(nx,ny,amass,pre1,pre1_in,vertx,verty,xc,yc, &
            bnd_type,Vertex_Type_norm_pre,mass_p,fx_pre,fy_pre)
   implicit none 
   integer   :: nx,ny
   real*8    :: amass(1:nx,1:ny),pre1(4,1:nx,1:ny),pre1_in(4,1:nx,1:ny)
   real*8    :: vertx(0:nx,0:ny),verty(0:nx,0:ny),mass_p(0:nx,0:ny)
   real*8    :: xc(1:nx,1:ny),yc(1:nx,1:ny)
   real*8    :: fx_pre(0:nx,0:ny),fy_pre(0:nx,0:ny)
	integer   :: bnd_type(0:nx,0:ny)
	real*8    :: Vertex_Type_norm_pre(0:nx,0:ny,3)   
   
   integer   :: i,j,k,k1,k2,k3,k4
   real*8    :: x1,x2,y1,y2,alength,anormalx,anormaly
   
   fx_pre = 0.d0 
   fy_pre = 0.d0    
   mass_p = 0.d0 

   ! --------------------------------
   ! |             |                |
   ! |             |                |
   ! |   *----3----*----2----*      |
   ! |   |         |         |      |
   ! |   4         |         1      |
   ! |   |         |         |      |
   ! ----*-------[i,j]-------*-------
   ! |   |         |         |      |
   ! |   5         |         8      |
   ! |   |         |         |      |
   ! |   *----6----*----7----*      |
   ! |             |                |
   ! |             |                |   
   ! --------------------------------
	do j = 1 , ny - 1
	do i = 1 , nx - 1
		! mass_p(i,j) = ( amass(i+1,j+1) + amass(i+1,j) + amass(i,j+1) + amass(i,j) ) / 4.d0
      
      do k = 1 , 8
         select case(k)
         case(1)
            x1 = (vertx(i,j)+vertx(i+1,j))/2.d0 ; y1 = (verty(i,j)+verty(i+1,j))/2.d0
            x2 = xc(i+1,j+1) ;	y2 = yc(i+1,j+1)
            k1 = 1 ; k2 = i+1 ; k3 = j+1 ; k4 = 1
         case(2)
            x1 = xc(i+1,j+1) ; y1 = yc(i+1,j+1)
            x2 = (vertx(i,j)+vertx(i,j+1))/2.d0 ; y2 = (verty(i,j)+verty(i,j+1))/2.d0
            k1 = 1 ; k2 = i+1 ; k3 = j+1 ; k4 = 4
         case(3)
            x1 = (vertx(i,j)+vertx(i,j+1))/2.d0 ;y1 = (verty(i,j)+verty(i,j+1))/2.d0
            x2 = xc(i,j+1);y2 = yc(i,j+1)
            k1 = 2 ; k2 = i ; k3 = j+1 ; k4 = 2
         case(4)
            x1 = xc(i,j+1) ;y1 = yc(i,j+1)
            x2 = (vertx(i,j)+vertx(i-1,j))/2.d0 ; y2 = (verty(i,j)+verty(i-1,j))/2.d0
            k1 = 2 ; k2 = i ; k3 = j+1 ; k4 = 1
         case(5)
            x1 = (vertx(i,j)+vertx(i-1,j))/2.d0 ;y1 = (verty(i,j)+verty(i-1,j))/2.d0
            x2 = xc(i,j) ; y2 = yc(i,j)
            k1 = 3 ; k2 = i ; k3 = j ; k4 = 3
         case(6)
            x1 = xc(i,j) ; y1 = yc(i,j)
            x2 = (vertx(i,j)+vertx(i,j-1))/2.d0 ;y2 = (verty(i,j)+verty(i,j-1))/2.d0
            k1 = 3 ; k2 = i ; k3 = j ; k4 = 2
         case(7)
            x1 = (vertx(i,j)+vertx(i,j-1))/2.d0 ;y1 = (verty(i,j)+verty(i,j-1))/2.d0
            x2 = xc(i+1,j) ; y2 = yc(i+1,j)
            k1 = 4 ; k2 = i+1 ; k3 = j ; k4 = 4
         case(8)
            x1 = xc(i+1,j) ; y1 = yc(i+1,j)
            x2 = (vertx(i,j)+vertx(i+1,j))/2.d0 ;y2 = (verty(i,j)+verty(i+1,j))/2.d0
            k1 = 4 ; k2 = i+1 ; k3 = j ; k4 = 3
         case default
            stop 'error,rdafe'
         end select
         call cal_length_normal(x1,x2,y1,y2,alength,anormalx,anormaly)
         fx_pre(i,j) = fx_pre(i,j) + (pre1(k1,k2,k3)+pre1_in(k4,k2,k3))*alength*anormalx 
         fy_pre(i,j) = fy_pre(i,j) + (pre1(k1,k2,k3)+pre1_in(k4,k2,k3))*alength*anormaly 
      enddo
	enddo
	enddo


   ! --------------------------------
   ! |             |                |
   ! |             |                |
   ! |   *----3----*----2----*      |
   ! |   |         |         |      |
   ! |   4         |         1      |
   ! |   |         |         |      |
   ! ----*----b2-[i,j]--b1---*-------
	do j = 0 , 0
	do i = 1 , nx - 1
		! mass_p(i,j) = ( amass(i+1,j+1) + amass(i,j+1) ) / 4.d0
      do k = 1 , 4
         select case(k)
         case(1)
            x1 = (vertx(i,j)+vertx(i+1,j))/2.d0 ; y1 = (verty(i,j)+verty(i+1,j))/2.d0
            x2 = xc(i+1,j+1) ;	y2 = yc(i+1,j+1)
            k1 = 1 ; k2 = i+1 ; k3 = j+1 ; k4 = 1
         case(2)
            x1 = xc(i+1,j+1) ; y1 = yc(i+1,j+1)
            x2 = (vertx(i,j)+vertx(i,j+1))/2.d0 ; y2 = (verty(i,j)+verty(i,j+1))/2.d0
            k1 = 1 ; k2 = i+1 ; k3 = j+1 ; k4 = 4
         case(3)
            x1 = (vertx(i,j)+vertx(i,j+1))/2.d0 ;y1 = (verty(i,j)+verty(i,j+1))/2.d0
            x2 = xc(i,j+1);y2 = yc(i,j+1)
            k1 = 2 ; k2 = i ; k3 = j+1 ; k4 = 2
         case(4)
            x1 = xc(i,j+1) ;y1 = yc(i,j+1)
            x2 = (vertx(i,j)+vertx(i-1,j))/2.d0 ; y2 = (verty(i,j)+verty(i-1,j))/2.d0
            k1 = 2 ; k2 = i ; k3 = j+1 ; k4 = 1
         case default
            stop 'error,rdafe'
         end select
         call cal_length_normal(x1,x2,y1,y2,alength,anormalx,anormaly)
         fx_pre(i,j) = fx_pre(i,j) + (pre1(k1,k2,k3)+pre1_in(k4,k2,k3))*alength*anormalx 
         fy_pre(i,j) = fy_pre(i,j) + (pre1(k1,k2,k3)+pre1_in(k4,k2,k3))*alength*anormaly 
      enddo      
      
		!> b1
		x1 = vertx(i,j) ;	y1 = verty(i,j)
		x2 = (vertx(i,j)+vertx(i+1,j))/2.d0 ; y2 = (verty(i,j)+verty(i+1,j))/2.d0
		call cal_length_normal(x1,x2,y1,y2,alength,anormalx,anormaly)
		fx_pre(i,j) = fx_pre(i,j) + Vertex_Type_norm_pre(i,j,3)*alength*anormalx 
		fy_pre(i,j) = fy_pre(i,j) + Vertex_Type_norm_pre(i,j,3)*alength*anormaly 

		!> b2
		x1 = (vertx(i,j)+vertx(i-1,j))/2.d0 ; y1 = (verty(i,j)+verty(i-1,j))/2.d0
		x2 = vertx(i,j) ;	y2 = verty(i,j)
		call cal_length_normal(x1,x2,y1,y2,alength,anormalx,anormaly)
		fx_pre(i,j) = fx_pre(i,j) + Vertex_Type_norm_pre(i,j,3)*alength*anormalx  
		fy_pre(i,j) = fy_pre(i,j) + Vertex_Type_norm_pre(i,j,3)*alength*anormaly     

	enddo
	enddo

   ! ----*----b2-[i,j]--b1---*-------
   ! |   |         |         |      |
   ! |   5         |         8      |
   ! |   |         |         |      |
   ! |   *----6----*----7----*      |
   ! |             |                |
   ! |             |                |   
   ! --------------------------------   
	do j = ny , ny
	do i = 1 , nx - 1
		! mass_p(i,j) = ( amass(i,j) + amass(i+1,j) ) / 4.d0
      do k = 5 , 8
         select case(k)
         case(5)
            x1 = (vertx(i,j)+vertx(i-1,j))/2.d0 ;y1 = (verty(i,j)+verty(i-1,j))/2.d0
            x2 = xc(i,j) ; y2 = yc(i,j)
            k1 = 3 ; k2 = i ; k3 = j ; k4 = 3
         case(6)
            x1 = xc(i,j) ; y1 = yc(i,j)
            x2 = (vertx(i,j)+vertx(i,j-1))/2.d0 ;y2 = (verty(i,j)+verty(i,j-1))/2.d0
            k1 = 3 ; k2 = i ; k3 = j ; k4 = 2
         case(7)
            x1 = (vertx(i,j)+vertx(i,j-1))/2.d0 ;y1 = (verty(i,j)+verty(i,j-1))/2.d0
            x2 = xc(i+1,j) ; y2 = yc(i+1,j)
            k1 = 4 ; k2 = i+1 ; k3 = j ; k4 = 4
         case(8)
            x1 = xc(i+1,j) ; y1 = yc(i+1,j)
            x2 = (vertx(i,j)+vertx(i+1,j))/2.d0 ;y2 = (verty(i,j)+verty(i+1,j))/2.d0
            k1 = 4 ; k2 = i+1 ; k3 = j ; k4 = 3
         case default
            stop 'error,rdafe'
         end select
         call cal_length_normal(x1,x2,y1,y2,alength,anormalx,anormaly)
         fx_pre(i,j) = fx_pre(i,j) + (pre1(k1,k2,k3)+pre1_in(k4,k2,k3))*alength*anormalx 
         fy_pre(i,j) = fy_pre(i,j) + (pre1(k1,k2,k3)+pre1_in(k4,k2,k3))*alength*anormaly 
      enddo     

		!> b1
		x1 = (vertx(i,j)+vertx(i+1,j))/2.d0 ; y1 = (verty(i,j)+verty(i+1,j))/2.d0
		x2 = vertx(i,j) ;	y2 = verty(i,j)
		call cal_length_normal(x1,x2,y1,y2,alength,anormalx,anormaly)
		fx_pre(i,j) = fx_pre(i,j) + Vertex_Type_norm_pre(i,j,3)*alength*anormalx 
		fy_pre(i,j) = fy_pre(i,j) + Vertex_Type_norm_pre(i,j,3)*alength*anormaly 

		!> b2
		x1 = vertx(i,j) ;	y1 = verty(i,j)
		x2 = (vertx(i,j)+vertx(i-1,j))/2.d0 ; y2 = (verty(i,j)+verty(i-1,j))/2.d0
		call cal_length_normal(x1,x2,y1,y2,alength,anormalx,anormaly)
		fx_pre(i,j) = fx_pre(i,j) + Vertex_Type_norm_pre(i,j,3)*alength*anormalx  
		fy_pre(i,j) = fy_pre(i,j) + Vertex_Type_norm_pre(i,j,3)*alength*anormaly           
       
	enddo
	enddo

   !   ------------------
   !   |                |
   !   |                |
   !   *----2----*      |
   !   |         |      |
   !   b1        1      |
   !   |         |      |
   ! [i,j]-------*-------
   !   |         |      |
   !   b2        8      |
   !   |         |      |
   !   *----7----*      |
   !   |                |
   !   |                |   
   !   ------------------   
	do j = 1 , ny - 1
	do i =  0 , 0
		! mass_p(i,j) = ( amass(i+1,j+1) + amass(i+1,j) ) / 4.d0
      do k = 1 , 4
         select case(k)
         case(1)
            x1 = (vertx(i,j)+vertx(i+1,j))/2.d0 ; y1 = (verty(i,j)+verty(i+1,j))/2.d0
            x2 = xc(i+1,j+1) ;	y2 = yc(i+1,j+1)
            k1 = 1 ; k2 = i+1 ; k3 = j+1 ; k4 = 1
         case(2)
            x1 = xc(i+1,j+1) ; y1 = yc(i+1,j+1)
            x2 = (vertx(i,j)+vertx(i,j+1))/2.d0 ; y2 = (verty(i,j)+verty(i,j+1))/2.d0
            k1 = 1 ; k2 = i+1 ; k3 = j+1 ; k4 = 4
         case(3)
            x1 = (vertx(i,j)+vertx(i,j-1))/2.d0 ; y1 = (verty(i,j)+verty(i,j-1))/2.d0
            x2 = xc(i+1,j) ; y2 = yc(i+1,j)
            k1 = 4 ; k2 = i+1 ; k3 = j ; k4 = 4
         case(4)
            x1 = xc(i+1,j) ; y1 = yc(i+1,j)
            x2 = (vertx(i,j)+vertx(i+1,j))/2.d0 ; y2 = (verty(i,j)+verty(i+1,j))/2.d0
            k1 = 4 ; k2 = i+1 ; k3 = j ; k4 = 3
         case default
            stop 'error,rdafe'
         end select
         call cal_length_normal(x1,x2,y1,y2,alength,anormalx,anormaly)
         fx_pre(i,j) = fx_pre(i,j) + (pre1(k1,k2,k3)+pre1_in(k4,k2,k3))*alength*anormalx 
         fy_pre(i,j) = fy_pre(i,j) + (pre1(k1,k2,k3)+pre1_in(k4,k2,k3))*alength*anormaly 
      !print*, i,j,fx_pre(i,j),fy_pre(i,j) ; pause 232
      enddo      
		!> b1
		x1 = (vertx(i,j)+vertx(i,j+1))/2.d0 ; y1 = (verty(i,j)+verty(i,j+1))/2.d0
		x2 = vertx(i,j) ;	y2 = verty(i,j)
		call cal_length_normal(x1,x2,y1,y2,alength,anormalx,anormaly)
		fx_pre(i,j) = fx_pre(i,j) + Vertex_Type_norm_pre(i,j,3)*alength*anormalx 
		fy_pre(i,j) = fy_pre(i,j) + Vertex_Type_norm_pre(i,j,3)*alength*anormaly 

		!> b2
		x1 = vertx(i,j) ;	y1 = verty(i,j)
		x2 = (vertx(i,j)+vertx(i,j-1))/2.d0 ; y2 = (verty(i,j)+verty(i,j-1))/2.d0
		call cal_length_normal(x1,x2,y1,y2,alength,anormalx,anormaly)
		fx_pre(i,j) = fx_pre(i,j) + Vertex_Type_norm_pre(i,j,3)*alength*anormalx  
		fy_pre(i,j) = fy_pre(i,j) + Vertex_Type_norm_pre(i,j,3)*alength*anormaly 
      
	enddo
	enddo

   ! ---------------
   ! |             |
   ! |             |
   ! |   *----3----*
   ! |   |         |
   ! |   4         b1
   ! |   |         |
   ! ----*-------[i,j]
   ! |   |         |
   ! |   5         b2
   ! |   |         |
   ! |   *----6----*
   ! |             |
   ! |             | 
   ! --------------- 
	do j = 1 , ny - 1
	do i =  nx , nx 
		! mass_p(i,j) = ( amass(i,j+1) + amass(i,j) ) / 4.d0
      do k = 3 , 6
         select case(k)
         case(3)
            x1 = (vertx(i,j)+vertx(i,j+1))/2.d0 ;y1 = (verty(i,j)+verty(i,j+1))/2.d0
            x2 = xc(i,j+1);y2 = yc(i,j+1)
            k1 = 2 ; k2 = i ; k3 = j+1 ; k4 = 2
         case(4)
            x1 = xc(i,j+1) ;y1 = yc(i,j+1)
            x2 = (vertx(i,j)+vertx(i-1,j))/2.d0 ; y2 = (verty(i,j)+verty(i-1,j))/2.d0
            k1 = 2 ; k2 = i ; k3 = j+1 ; k4 = 1
         case(5)
            x1 = (vertx(i,j)+vertx(i-1,j))/2.d0 ;y1 = (verty(i,j)+verty(i-1,j))/2.d0
            x2 = xc(i,j) ; y2 = yc(i,j)
            k1 = 3 ; k2 = i ; k3 = j ; k4 = 3
         case(6)
            x1 = xc(i,j) ; y1 = yc(i,j)
            x2 = (vertx(i,j)+vertx(i,j-1))/2.d0 ;y2 = (verty(i,j)+verty(i,j-1))/2.d0
            k1 = 3 ; k2 = i ; k3 = j ; k4 = 2
         case default
            stop 'error,rdafe'
         end select
         call cal_length_normal(x1,x2,y1,y2,alength,anormalx,anormaly)
         fx_pre(i,j) = fx_pre(i,j) + (pre1(k1,k2,k3)+pre1_in(k4,k2,k3))*alength*anormalx 
         fy_pre(i,j) = fy_pre(i,j) + (pre1(k1,k2,k3)+pre1_in(k4,k2,k3))*alength*anormaly 
      enddo
		!> b1
		x1 = vertx(i,j) ;	y1 = verty(i,j)
		x2 = (vertx(i,j)+vertx(i,j+1))/2.d0 ; y2 = (verty(i,j)+verty(i,j+1))/2.d0
		call cal_length_normal(x1,x2,y1,y2,alength,anormalx,anormaly)
		fx_pre(i,j) = fx_pre(i,j) + Vertex_Type_norm_pre(i,j,3)*alength*anormalx 
		fy_pre(i,j) = fy_pre(i,j) + Vertex_Type_norm_pre(i,j,3)*alength*anormaly 
		!> b2
		x1 = (vertx(i,j)+vertx(i,j-1))/2.d0 ; y1 = (verty(i,j)+verty(i,j-1))/2.d0
		x2 = vertx(i,j) ;	y2 = verty(i,j)
		call cal_length_normal(x1,x2,y1,y2,alength,anormalx,anormaly)
		fx_pre(i,j) = fx_pre(i,j) + Vertex_Type_norm_pre(i,j,3)*alength*anormalx  
		fy_pre(i,j) = fy_pre(i,j) + Vertex_Type_norm_pre(i,j,3)*alength*anormaly

	enddo
	enddo

   !   ------------------
   !   |                |
   !   |                |
   !   *----2----*      |
   !   |         |      |
   !   b2        1      |
   !   |         |      |
   ! [i,j]--b1---*-------
 
	do j = 0 , 0
	do i = 0 , 0
		! mass_p(i,j) = ( amass(i+1,j+1) ) / 4.d0
      do k = 1 , 2
         select case(k)
         case(1)
            x1 = (vertx(i,j)+vertx(i+1,j))/2.d0 ; y1 = (verty(i,j)+verty(i+1,j))/2.d0
            x2 = xc(i+1,j+1) ;	y2 = yc(i+1,j+1)
            k1 = 1 ; k2 = i+1 ; k3 = j+1 ; k4 = 1
         case(2)
            x1 = xc(i+1,j+1) ; y1 = yc(i+1,j+1)
            x2 = (vertx(i,j)+vertx(i,j+1))/2.d0 ; y2 = (verty(i,j)+verty(i,j+1))/2.d0
            k1 = 1 ; k2 = i+1 ; k3 = j+1 ; k4 = 4
         case default
            stop 'error,rdafe'
         end select
         call cal_length_normal(x1,x2,y1,y2,alength,anormalx,anormaly)
         fx_pre(i,j) = fx_pre(i,j) + (pre1(k1,k2,k3)+pre1_in(k4,k2,k3))*alength*anormalx 
         fy_pre(i,j) = fy_pre(i,j) + (pre1(k1,k2,k3)+pre1_in(k4,k2,k3))*alength*anormaly 
      enddo
      
      if ( bnd_type(i,j) == 241 ) then 
         !> b1
         x1 = vertx(i,j) ;	y1 = verty(i,j)
         x2 = (vertx(i,j)+vertx(i+1,j))/2.d0 ; y2 = (verty(i,j)+verty(i+1,j))/2.d0
         call cal_length_normal(x1,x2,y1,y2,alength,anormalx,anormaly)
         fx_pre(i,j) = fx_pre(i,j) + Vertex_Type_norm_pre(i,j,3)*alength*anormalx 
         fy_pre(i,j) = fy_pre(i,j) + Vertex_Type_norm_pre(i,j,3)*alength*anormaly 

      else if ( bnd_type(i,j) == 242 ) then 
         !> b2
         x1 = (vertx(i,j)+vertx(i,j+1))/2.d0 ; y1 = (verty(i,j)+verty(i,j+1))/2.d0
         x2 = vertx(i,j) ;	y2 = verty(i,j)
         call cal_length_normal(x1,x2,y1,y2,alength,anormalx,anormaly)
         fx_pre(i,j) = fx_pre(i,j) + Vertex_Type_norm_pre(i,j,3)*alength*anormalx  
         fy_pre(i,j) = fy_pre(i,j) + Vertex_Type_norm_pre(i,j,3)*alength*anormaly       
         
      else 
         !> b1
         x1 = vertx(i,j) ;	y1 = verty(i,j)
         x2 = (vertx(i,j)+vertx(i+1,j))/2.d0 ; y2 = (verty(i,j)+verty(i+1,j))/2.d0
         call cal_length_normal(x1,x2,y1,y2,alength,anormalx,anormaly)
         fx_pre(i,j) = fx_pre(i,j) + Vertex_Type_norm_pre(i,j,3)*alength*anormalx 
         fy_pre(i,j) = fy_pre(i,j) + Vertex_Type_norm_pre(i,j,3)*alength*anormaly       
      
         !> b2
         x1 = (vertx(i,j)+vertx(i,j+1))/2.d0 ; y1 = (verty(i,j)+verty(i,j+1))/2.d0
         x2 = vertx(i,j) ;	y2 = verty(i,j)
         call cal_length_normal(x1,x2,y1,y2,alength,anormalx,anormaly)
         fx_pre(i,j) = fx_pre(i,j) + Vertex_Type_norm_pre(i,j,3)*alength*anormalx  
         fy_pre(i,j) = fy_pre(i,j) + Vertex_Type_norm_pre(i,j,3)*alength*anormaly
         
      endif 
      
            
	enddo
	enddo

   ! ---------------
   ! |             |  
   ! |             |  
   ! |   *----3----*
   ! |   |         |  
   ! |   4         b1  
   ! |   |         |  
   ! ----*---b2--[i,j]
	do j = 0 , 0
	do i = nx , nx 
		! mass_p(i,j) = ( amass(i,j+1) ) / 4.d0
      do k = 3 , 4
         select case(k)
         case(3)
            x1 = (vertx(i,j)+vertx(i,j+1))/2.d0 ;y1 = (verty(i,j)+verty(i,j+1))/2.d0
            x2 = xc(i,j+1);y2 = yc(i,j+1)
            k1 = 2 ; k2 = i ; k3 = j+1 ; k4 = 2
         case(4)
            x1 = xc(i,j+1) ;y1 = yc(i,j+1)
            x2 = (vertx(i,j)+vertx(i-1,j))/2.d0 ; y2 = (verty(i,j)+verty(i-1,j))/2.d0
            k1 = 2 ; k2 = i ; k3 = j+1 ; k4 = 1
         case default
            stop 'error,rdafe'
         end select
         call cal_length_normal(x1,x2,y1,y2,alength,anormalx,anormaly)
         fx_pre(i,j) = fx_pre(i,j) + (pre1(k1,k2,k3)+pre1_in(k4,k2,k3))*alength*anormalx 
         fy_pre(i,j) = fy_pre(i,j) + (pre1(k1,k2,k3)+pre1_in(k4,k2,k3))*alength*anormaly 
      enddo
      
      if ( bnd_type(i,j) == 241 ) then       
         !> b1 
         x1 = vertx(i,j) ;	y1 = verty(i,j)
         x2 = (vertx(i,j)+vertx(i,j+1))/2.d0 ; y2 = (verty(i,j)+verty(i,j+1))/2.d0
         call cal_length_normal(x1,x2,y1,y2,alength,anormalx,anormaly)
         fx_pre(i,j) = fx_pre(i,j) + Vertex_Type_norm_pre(i,j,3)*alength*anormalx 
         fy_pre(i,j) = fy_pre(i,j) + Vertex_Type_norm_pre(i,j,3)*alength*anormaly       
      
      else if ( bnd_type(i,j) == 242 ) then    
         !> b2 
         x1 = (vertx(i,j)+vertx(i-1,j))/2.d0 ; y1 = (verty(i,j)+verty(i-1,j))/2.d0
         x2 = vertx(i,j) ;	y2 = verty(i,j)
         call cal_length_normal(x1,x2,y1,y2,alength,anormalx,anormaly)
         fx_pre(i,j) = fx_pre(i,j) + Vertex_Type_norm_pre(i,j,3)*alength*anormalx  
         fy_pre(i,j) = fy_pre(i,j) + Vertex_Type_norm_pre(i,j,3)*alength*anormaly        
      else  
         !> b1 
         x1 = vertx(i,j) ;	y1 = verty(i,j)
         x2 = (vertx(i,j)+vertx(i,j+1))/2.d0 ; y2 = (verty(i,j)+verty(i,j+1))/2.d0
         call cal_length_normal(x1,x2,y1,y2,alength,anormalx,anormaly)
         fx_pre(i,j) = fx_pre(i,j) + Vertex_Type_norm_pre(i,j,3)*alength*anormalx 
         fy_pre(i,j) = fy_pre(i,j) + Vertex_Type_norm_pre(i,j,3)*alength*anormaly 
         !> b2 
         x1 = (vertx(i,j)+vertx(i-1,j))/2.d0 ; y1 = (verty(i,j)+verty(i-1,j))/2.d0
         x2 = vertx(i,j) ;	y2 = verty(i,j)
         call cal_length_normal(x1,x2,y1,y2,alength,anormalx,anormaly)
         fx_pre(i,j) = fx_pre(i,j) + Vertex_Type_norm_pre(i,j,3)*alength*anormalx  
         fy_pre(i,j) = fy_pre(i,j) + Vertex_Type_norm_pre(i,j,3)*alength*anormaly            
      
      endif

	enddo
	enddo	


   ! [i,j]--b2---*-------
   !   |         |      |
   !   b1        8      |
   !   |         |      |
   !   *----7----*      |
   !   |                |
   !   |                |   
   !   ------------------   
	do j = ny , ny
	do i = 0 , 0
		! mass_p(i,j) = ( amass(i+1,j) ) / 4.d0
      do k = 7 , 8
         select case(k)
         case(7)
            x1 = (vertx(i,j)+vertx(i,j-1))/2.d0 ;y1 = (verty(i,j)+verty(i,j-1))/2.d0
            x2 = xc(i+1,j) ; y2 = yc(i+1,j)
            k1 = 4 ; k2 = i+1 ; k3 = j ; k4 = 4
         case(8)
            x1 = xc(i+1,j) ; y1 = yc(i+1,j)
            x2 = (vertx(i,j)+vertx(i+1,j))/2.d0 ;y2 = (verty(i,j)+verty(i+1,j))/2.d0
            k1 = 4 ; k2 = i+1 ; k3 = j ; k4 = 3
         case default
            stop 'error,rdafe'
         end select
         call cal_length_normal(x1,x2,y1,y2,alength,anormalx,anormaly)
         fx_pre(i,j) = fx_pre(i,j) + (pre1(k1,k2,k3)+pre1_in(k4,k2,k3))*alength*anormalx 
         fy_pre(i,j) = fy_pre(i,j) + (pre1(k1,k2,k3)+pre1_in(k4,k2,k3))*alength*anormaly 
      enddo
      
      if ( bnd_type(i,j) == 241 ) then           
         !> b1
         x1 = vertx(i,j) ;	y1 = verty(i,j)
         x2 = (vertx(i,j)+vertx(i,j-1))/2.d0 ; y2 = (verty(i,j)+verty(i,j-1))/2.d0
         call cal_length_normal(x1,x2,y1,y2,alength,anormalx,anormaly)
         fx_pre(i,j) = fx_pre(i,j) + Vertex_Type_norm_pre(i,j,3)*alength*anormalx  
         fy_pre(i,j) = fy_pre(i,j) + Vertex_Type_norm_pre(i,j,3)*alength*anormaly      
      else if ( bnd_type(i,j) == 242 ) then     
         !> b2
         x1 = (vertx(i,j)+vertx(i+1,j))/2.d0 ; y1 = (verty(i,j)+verty(i+1,j))/2.d0
         x2 = vertx(i,j) ;	y2 = verty(i,j)
         call cal_length_normal(x1,x2,y1,y2,alength,anormalx,anormaly)
         fx_pre(i,j) = fx_pre(i,j) + Vertex_Type_norm_pre(i,j,3)*alength*anormalx 
         fy_pre(i,j) = fy_pre(i,j) + Vertex_Type_norm_pre(i,j,3)*alength*anormaly 
      else
         !> b1
         x1 = vertx(i,j) ;	y1 = verty(i,j)
         x2 = (vertx(i,j)+vertx(i,j-1))/2.d0 ; y2 = (verty(i,j)+verty(i,j-1))/2.d0
         call cal_length_normal(x1,x2,y1,y2,alength,anormalx,anormaly)
         fx_pre(i,j) = fx_pre(i,j) + Vertex_Type_norm_pre(i,j,3)*alength*anormalx  
         fy_pre(i,j) = fy_pre(i,j) + Vertex_Type_norm_pre(i,j,3)*alength*anormaly        
         !> b2
         x1 = (vertx(i,j)+vertx(i+1,j))/2.d0 ; y1 = (verty(i,j)+verty(i+1,j))/2.d0
         x2 = vertx(i,j) ;	y2 = verty(i,j)
         call cal_length_normal(x1,x2,y1,y2,alength,anormalx,anormaly)
         fx_pre(i,j) = fx_pre(i,j) + Vertex_Type_norm_pre(i,j,3)*alength*anormalx 
         fy_pre(i,j) = fy_pre(i,j) + Vertex_Type_norm_pre(i,j,3)*alength*anormaly 
      endif 

	enddo
	enddo
   
   ! ----*----b1-[i,j]
   ! |   |         |  
   ! |   5         b2  
   ! |   |         |  
   ! |   *----6----*
   ! |             |  
   ! |             |   
   ! ---------------    
   
	do j = ny , ny
	do i = nx , nx
		! mass_p(i,j) = ( amass(i,j) ) / 4.d0
      do k = 5 , 6
         select case(k)
         case(5)
            x1 = (vertx(i,j)+vertx(i-1,j))/2.d0 ;y1 = (verty(i,j)+verty(i-1,j))/2.d0
            x2 = xc(i,j) ; y2 = yc(i,j)
            k1 = 3 ; k2 = i ; k3 = j ; k4 = 3
         case(6)
            x1 = xc(i,j) ; y1 = yc(i,j)
            x2 = (vertx(i,j)+vertx(i,j-1))/2.d0 ;y2 = (verty(i,j)+verty(i,j-1))/2.d0
            k1 = 3 ; k2 = i ; k3 = j ; k4 = 2
         case default
            stop 'error,rdafe'
         end select
         call cal_length_normal(x1,x2,y1,y2,alength,anormalx,anormaly)
         fx_pre(i,j) = fx_pre(i,j) + (pre1(k1,k2,k3)+pre1_in(k4,k2,k3))*alength*anormalx 
         fy_pre(i,j) = fy_pre(i,j) + (pre1(k1,k2,k3)+pre1_in(k4,k2,k3))*alength*anormaly 
      enddo

      if ( bnd_type(i,j) == 241 ) then           
         !> b1
         x1 = vertx(i,j) ;	y1 = verty(i,j)
         x2 = (vertx(i,j)+vertx(i-1,j))/2.d0 ; y2 = (verty(i,j)+verty(i-1,j))/2.d0
         call cal_length_normal(x1,x2,y1,y2,alength,anormalx,anormaly)
         fx_pre(i,j) = fx_pre(i,j) + Vertex_Type_norm_pre(i,j,3)*alength*anormalx  
         fy_pre(i,j) = fy_pre(i,j) + Vertex_Type_norm_pre(i,j,3)*alength*anormaly  
      else if ( bnd_type(i,j) == 242 ) then     
         !> b2
         x1 = (vertx(i,j)+vertx(i,j-1))/2.d0 ; y1 = (verty(i,j)+verty(i,j-1))/2.d0
         x2 = vertx(i,j) ;	y2 = verty(i,j)
         call cal_length_normal(x1,x2,y1,y2,alength,anormalx,anormaly)
         fx_pre(i,j) = fx_pre(i,j) + Vertex_Type_norm_pre(i,j,3)*alength*anormalx 
         fy_pre(i,j) = fy_pre(i,j) + Vertex_Type_norm_pre(i,j,3)*alength*anormaly 
      else 
         !> b1
         x1 = vertx(i,j) ;	y1 = verty(i,j)
         x2 = (vertx(i,j)+vertx(i-1,j))/2.d0 ; y2 = (verty(i,j)+verty(i-1,j))/2.d0
         call cal_length_normal(x1,x2,y1,y2,alength,anormalx,anormaly)
         fx_pre(i,j) = fx_pre(i,j) + Vertex_Type_norm_pre(i,j,3)*alength*anormalx  
         fy_pre(i,j) = fy_pre(i,j) + Vertex_Type_norm_pre(i,j,3)*alength*anormaly        
         !> b2
         x1 = (vertx(i,j)+vertx(i,j-1))/2.d0 ; y1 = (verty(i,j)+verty(i,j-1))/2.d0
         x2 = vertx(i,j) ;	y2 = verty(i,j)
         call cal_length_normal(x1,x2,y1,y2,alength,anormalx,anormaly)
         fx_pre(i,j) = fx_pre(i,j) + Vertex_Type_norm_pre(i,j,3)*alength*anormalx 
         fy_pre(i,j) = fy_pre(i,j) + Vertex_Type_norm_pre(i,j,3)*alength*anormaly 
      endif 
           
	enddo
   enddo   
   
end subroutine cal_node_force_cy
 
