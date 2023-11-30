


subroutine cal_mjsls_details(nx,ny,den,pre,sound,enin,eos,vertx,verty,amass, &
            gamar,vstarx,vstary,dt,bnd_type,Vertex_Type_norm_pre,g)
   implicit none 
   
	integer,intent(in) :: nx,ny
	real*8 :: den(1:nx,1:ny),pre(1:nx,1:ny),amass(1:nx,1:ny)
	real*8 :: sound(1:nx,1:ny),enin(1:nx,1:ny)
	real*8 :: xc(1:nx,1:ny),yc(1:nx,1:ny),g
   real*8 :: pre1(4,1:nx,1:ny),pre1_in(4,1:nx,1:ny)
	integer :: eos(1:nx,1:ny),bnd_type(0:nx,0:ny)
	real*8  :: Vertex_Type_norm_pre(0:nx,0:ny,3)
   
	real*8 :: vertx(0:nx,0:ny),verty(0:nx,0:ny)
	real*8 :: vstarx(0:nx,0:ny),vstary(0:nx,0:ny)
	real*8 :: velx(0:nx,0:ny),vely(0:nx,0:ny),mass_p(0:nx,0:ny)
   real*8 :: fx_pre(0:nx,0:ny),fy_pre(0:nx,0:ny)
	real*8  :: gamar,dt
   
	real*8 :: den_new(1:nx,1:ny),pre_new(1:nx,1:ny)
   real*8 :: sound_new(1:nx,1:ny),enin_new(1:nx,1:ny)
	real*8 :: vertx_new(0:nx,0:ny),verty_new(0:nx,0:ny)
	real*8 :: vstarx_new(0:nx,0:ny),vstary_new(0:nx,0:ny)

   real*8 :: pre2(4,1:nx,1:ny),pre2_in(4,1:nx,1:ny)
	real*8 :: xc_pre(1:nx,1:ny),yc_pre(1:nx,1:ny)
   real*8 :: fx_cor(0:nx,0:ny),fy_cor(0:nx,0:ny)

   call cal_mjsls_details_subcell_pre4_in(nx,ny,eos,vertx,verty,vstarx,vstary, &   
            den,pre,sound,gamar,amass,xc,yc,pre1,pre1_in)
            
   call cal_node_force(nx,ny,amass,pre1,pre1_in,vertx,verty,xc,yc, &
            bnd_type,Vertex_Type_norm_pre,mass_p,fx_pre,fy_pre)      

   call cal_node_velocity(nx,ny,mass_p,dt,bnd_type,Vertex_Type_norm_pre, &
         vstarx,vstary,vstarx_new,vstary_new,fx_pre,fy_pre,g)         

   velx = (vstarx + vstarx_new) / 2.d0 
   vely = (vstary + vstary_new) / 2.d0   
   call update_vertex(nx,ny,vertx,verty,vertx_new,verty_new,velx,vely,dt)

   call update_den_enin(nx,ny,eos,vertx,verty,vertx_new,verty_new,velx,vely,dt,&
            pre1,den,enin,amass,gamar,den_new,enin_new,pre_new,sound_new,g) 
    
   call cal_mjsls_details_subcell_pre4_in(nx,ny,eos,vertx_new,verty_new,vstarx_new,vstary_new, &   
            den_new,pre_new,sound_new,gamar,amass,xc_pre,yc_pre,pre2,pre2_in)

   call cal_node_force(nx,ny,amass,(pre1+pre2)/2.d0,(pre1_in+pre2_in)/2.d0,&
            vertx_new,verty_new,xc_pre,yc_pre, &
            bnd_type,Vertex_Type_norm_pre,mass_p,fx_cor,fy_cor)         

   call cal_node_velocity(nx,ny,mass_p,dt,bnd_type,Vertex_Type_norm_pre, &
         vstarx,vstary,vstarx_new,vstary_new,(fx_pre+fx_cor)/2.d0,(fy_pre+fy_cor)/2.d0,g)                

   velx = (vstarx + vstarx_new) / 2.d0 
   vely = (vstary + vstary_new) / 2.d0   
   call update_vertex(nx,ny,vertx,verty,vertx_new,verty_new,velx,vely,dt)

   call update_den_enin(nx,ny,eos,vertx,verty,vertx_new,verty_new,velx,vely,dt,&
            (pre1+pre2)/2.d0,den,enin,amass,gamar,den_new,enin_new,pre_new,sound_new,g)    
        
   vertx  = vertx_new
   verty  = verty_new 
   vstarx = vstarx_new
   vstary = vstary_new       
   den    = den_new
   pre    = pre_new
   sound  = sound_new
   enin   = enin_new    

end subroutine cal_mjsls_details

   
subroutine cal_node_force(nx,ny,amass,pre1,pre1_in,vertx,verty,xc,yc, &
            bnd_type,Vertex_Type_norm_pre,mass_p,fx_pre,fy_pre)
   implicit none 
   integer   :: nx,ny
   real*8    :: amass(1:nx,1:ny),pre1(4,1:nx,1:ny),pre1_in(4,1:nx,1:ny)
   real*8    :: vertx(0:nx,0:ny),verty(0:nx,0:ny)
   real*8    :: xc(1:nx,1:ny),yc(1:nx,1:ny)
   real*8    :: fx_pre(0:nx,0:ny),fy_pre(0:nx,0:ny),mass_p(0:nx,0:ny)
	integer   :: bnd_type(0:nx,0:ny)
	real*8    :: Vertex_Type_norm_pre(0:nx,0:ny,3)   
   
   integer   :: i,j,k,k1,k2,k3,k4
   real*8    :: x1,x2,y1,y2,alength,anormalx,anormaly
   
   fx_pre = 0.d0 
   fy_pre = 0.d0    

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
		mass_p(i,j) = ( amass(i+1,j+1) + amass(i+1,j) + amass(i,j+1) + amass(i,j) ) / 4.d0
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
		mass_p(i,j) = ( amass(i+1,j+1) + amass(i,j+1) ) / 4.d0
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
		mass_p(i,j) = ( amass(i,j) + amass(i+1,j) ) / 4.d0
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
		mass_p(i,j) = ( amass(i+1,j+1) + amass(i+1,j) ) / 4.d0
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
		mass_p(i,j) = ( amass(i,j+1) + amass(i,j) ) / 4.d0
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
		mass_p(i,j) = ( amass(i+1,j+1) ) / 4.d0
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
		mass_p(i,j) = ( amass(i,j+1) ) / 4.d0
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
		mass_p(i,j) = ( amass(i+1,j) ) / 4.d0
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
		mass_p(i,j) = ( amass(i,j) ) / 4.d0
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
   
end subroutine cal_node_force
 

   
subroutine cal_node_velocity(nx,ny,mass_p,dt,bnd_type,Vertex_Type_norm_pre, &
            vstarx,vstary,vstarx_pre,vstary_pre,fx_pre,fy_pre,g)
            
   implicit none 
   integer   :: nx,ny
   real*8    :: vstarx(0:nx,0:ny),vstary(0:nx,0:ny),mass_p(0:nx,0:ny)
   real*8    :: vstarx_pre(0:nx,0:ny),vstary_pre(0:nx,0:ny)
   real*8    :: fx_pre(0:nx,0:ny),fy_pre(0:nx,0:ny),dt,g
	integer   :: bnd_type(0:nx,0:ny)
	real*8    :: Vertex_Type_norm_pre(0:nx,0:ny,3)   
   
   integer  :: i,j
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
   do j = 0 , ny 
   do i = 0 , nx 
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
   enddo
   enddo 


   do j = 1 , ny-1 
   do i = 0 , nx    
      vstary_pre(i,j) = vstary_pre(i,j) - dt * g
   enddo
   enddo 

end subroutine cal_node_velocity
 

subroutine cal_mjsls_details_subcell_pre4_in(nx,ny,eos,vertx,verty,vstarx,vstary, &   
            den,pre,sound,gamar,amass,xc,yc,pre1,pre1_in)
	implicit none 
   integer :: nx,ny
	real*8 :: den(1:nx,1:ny),pre(1:nx,1:ny)
	real*8 :: sound(1:nx,1:ny),amass(1:nx,1:ny)
	real*8 :: vertx(0:nx,0:ny),verty(0:nx,0:ny)
	real*8 :: vstarx(0:nx,0:ny),vstary(0:nx,0:ny)
   integer :: eos(1:nx,1:ny)
	real*8  :: gamar 
	real*8,intent(out)  :: pre1(4,1:nx,1:ny),pre1_in(4,1:nx,1:ny)
	real*8,intent(out)  :: xc(1:nx,1:ny),yc(1:nx,1:ny)
   

   integer  :: i,j 
   real*8   :: phi(1:nx,1:ny)

   ! call cal_limiter_phi( nx,ny,vertx,verty,vstarx,vstary,phi )

   !> compute pre1 and pre1_in
	do j = 1 , ny 
	do i = 1 , nx 
      ! print*, i,j
      call cal_eos(eos(i,j),gamar)
		call cal_subcell_p4_in(                                            &
				vertx(i-1,j-1),verty(i-1,j-1),vertx(i,j-1),verty(i,j-1),     &
				vertx(i,j),verty(i,j),vertx(i-1,j),verty(i-1,j),             &
				vstarx(i-1,j-1),vstary(i-1,j-1),vstarx(i,j-1),vstary(i,j-1), &
				vstarx(i,j),vstary(i,j),vstarx(i-1,j),vstary(i-1,j),         &
            xc(i,j),yc(i,j),den(i,j),pre(i,j),sound(i,j),pre1(1:4,i,j),  &
            gamar,amass(i,j),pre1_in(1:4,i,j),phi(i,j) )            
	enddo
   enddo
   ! pause 798
	! do j = 1 , ny 
	! do i = 1 , 1    
   !    pre1(1:4,i,j) = pre(i,j)
   ! enddo
   ! enddo 

end subroutine cal_mjsls_details_subcell_pre4_in



subroutine cal_subcell_p4_in(x1,y1,x2,y2,x3,y3,x4,y4,u1,v1,u2,v2,u3,v3,u4,v4,xc,yc, &
			den,pre,sound,pre1,gamar,amass,pre1_in,phi)
	implicit none 
	real*8,intent(in)   :: x1,y1,x2,y2,x3,y3,x4,y4
	real*8,intent(in)   :: u1,v1,u2,v2,u3,v3,u4,v4
	real*8,intent(in)   :: den,pre,sound,gamar,amass,phi
	real*8,intent(out)  :: pre1(4),pre1_in(4),xc,yc

	integer  :: i,i1,i2,side,n
	real*8   :: x12,y12,x23,y23,x34,y34,x41,y41,sub_area(4) 
	real*8   :: xa(4),ya(4),ua(4),va(4),alength,anormalx,anormaly
	real*8   :: delta_v,clength,area,xx1,yy1,hbar,delta_tau,xx2,yy2
   real*8   :: ustar,vstar,hbar4(4),pstar4(4),rwai,rnei,uc,vc,tden 

	xc = (x1+x2+x3+x4) / 4.d0 
	yc = (y1+y2+y3+y4) / 4.d0 

   xa(1) = x1   ; ya(1) = y1
   xa(2) = x2   ; ya(2) = y2
   xa(3) = x3   ; ya(3) = y3 
   xa(4) = x4   ; ya(4) = y4 
   ua(1) = u1   ; va(1) = v1
   ua(2) = u2   ; va(2) = v2
   ua(3) = u3   ; va(3) = v3 
   ua(4) = u4   ; va(4) = v4 
   uc = sum(ua(:)) / 4.d0 
   vc = sum(va(:)) / 4.d0 

   delta_v = 0.d0 
   area    = 0.d0 
   clength = 0.d0 
   do i = 1 , 4
      if( i == 4 ) then 
         i1 = 4 
         i2 = 1 
      else 
         i1 = i 
         i2 = i + 1
      endif 
      xx1 = xa(i1) ; yy1 = ya(i1)
      xx2 = xa(i2) ; yy2 = ya(i2)
      call cal_length_normal(xx1,xx2,yy1,yy2,alength,anormalx,anormaly)
      delta_v = delta_v + alength*( (ua(i1)+ua(i2))/2.d0*anormalx + &
                                    (va(i1)+va(i2))/2.d0*anormaly     )	
      area = area + 0.5d0*(xx1*yy2-xx2*yy1)
      clength = clength + alength
   enddo 
   
   hbar = area/clength * 2.d0 * 2.d0 

   delta_tau = delta_v * hbar / area / den / sound

   if( delta_tau < 0.d0 ) then 
         pre1 = pre - (den*sound)**2*delta_tau  &
               + (gamar+1.d0)/2.d0*den**3*sound**2*delta_tau**2 
   else 
      pre1 = pre
   endif  

	x12 = (x1+x2)/2.d0
	y12 = (y1+y2)/2.d0 
	x23 = (x2+x3)/2.d0
	y23 = (y2+y3)/2.d0 
	x34 = (x3+x4)/2.d0 
	y34 = (y3+y4)/2.d0 
	x41 = (x4+x1)/2.d0 
	y41 = (y4+y1)/2.d0  
 
 	n = 1
	do side = 1 , 4
		select case (side)
		case(1)
			xa(1) = x1   ; ya(1) = y1
			xa(2) = x12  ; ya(2) = y12
			xa(3) = xc   ; ya(3) = yc 
			xa(4) = x41  ; ya(4) = y41 
		case(2)
			xa(1) = x2   ; ya(1) = y2
			xa(2) = x23  ; ya(2) = y23
			xa(3) = xc   ; ya(3) = yc 
			xa(4) = x12  ; ya(4) = y12 
		case(3)
			xa(1) = x3   ; ya(1) = y3
			xa(2) = x34  ; ya(2) = y34
			xa(3) = xc   ; ya(3) = yc 
			xa(4) = x23  ; ya(4) = y23 
		case(4)
			xa(1) = x4   ; ya(1) = y4
			xa(2) = x41  ; ya(2) = y41
			xa(3) = xc   ; ya(3) = yc 
			xa(4) = x34  ; ya(4) = y34 
		case default
			stop 'error,fdaoe34'
      end select 
      area = 0.d0 
		do i = 1 , 4
         i1 = i
			if( i == 4 ) then 
				i2 = 1 
			else 
				i2 = i + 1
			endif 
			xx1 = xa(i1) ; yy1 = ya(i1)
			xx2 = xa(i2) ; yy2 = ya(i2)
			area = area + 0.5d0*(xx1*yy2-xx2*yy1)
      enddo 
      sub_area(side) = area 
   enddo
   do side = 1 , 4
      i1 = side
      if(side == 4) then
         i2 = 1
      else 
         i2 = side + 1
      endif 
      pre1_in(side)  = sound**2 * (amass/2.d0/(sub_area(i1)+sub_area(i2)) - den) 
   enddo 

end subroutine cal_subcell_p4_in


