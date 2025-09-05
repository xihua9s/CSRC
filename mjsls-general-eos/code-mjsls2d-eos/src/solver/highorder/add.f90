

subroutine cal_limiter_phi( nx,ny,vertx,verty,vstarx,vstary,phi )   
	implicit none 
   integer :: nx,ny
	real*8 :: vertx(0:nx,0:ny),verty(0:nx,0:ny)
	real*8 :: vstarx(0:nx,0:ny),vstary(0:nx,0:ny)
	real*8,intent(out)  :: phi(1:nx,1:ny)

	real*8   :: x1,y1,x2,y2,x3,y3,x4,y4
   real*8   :: u1,v1,u2,v2,u3,v3,u4,v4
   real*8   :: phi1,phi2,phi3,phi4
   integer  :: i,j 

   phi = 1.d0 

   do j = 2 , ny-1
   do i = 2 , nx-1
      x1 = vertx(i-2,j-1) ; y1 = verty(i-2,j-1)
      x2 = vertx(i-1,j-1) ; y2 = verty(i-1,j-1)
      x3 = vertx(i-0,j-1) ; y3 = verty(i-0,j-1)
      x4 = vertx(i+1,j-1) ; y4 = verty(i+1,j-1)
      u1 = vstarx(i-2,j-1) ; v1 = vstary(i-2,j-1)
      u2 = vstarx(i-1,j-1) ; v2 = vstary(i-1,j-1)
      u3 = vstarx(i-0,j-1) ; v3 = vstary(i-0,j-1)
      u4 = vstarx(i+1,j-1) ; v4 = vstary(i+1,j-1)
      call limiter_phi_4(x1,y1,x2,y2,x3,y3,x4,y4,&
            u1,v1,u2,v2,u3,v3,u4,v4,phi1)

      x1 = vertx(i-2,j-0) ; y1 = verty(i-2,j-0)
      x2 = vertx(i-1,j-0) ; y2 = verty(i-1,j-0)
      x3 = vertx(i-0,j-0) ; y3 = verty(i-0,j-0)
      x4 = vertx(i+1,j-0) ; y4 = verty(i+1,j-0)
      u1 = vstarx(i-2,j-0) ; v1 = vstary(i-2,j-0)
      u2 = vstarx(i-1,j-0) ; v2 = vstary(i-1,j-0)
      u3 = vstarx(i-0,j-0) ; v3 = vstary(i-0,j-0)
      u4 = vstarx(i+1,j-0) ; v4 = vstary(i+1,j-0)
      call limiter_phi_4(x1,y1,x2,y2,x3,y3,x4,y4,&
            u1,v1,u2,v2,u3,v3,u4,v4,phi2)

      x1 = vertx(i-1,j-2) ; y1 = verty(i-1,j-2)
      x2 = vertx(i-1,j-1) ; y2 = verty(i-1,j-1)
      x3 = vertx(i-1,j-0) ; y3 = verty(i-1,j-0)
      x4 = vertx(i-1,j+1) ; y4 = verty(i-1,j+1)
      u1 = vstarx(i-1,j-2) ; v1 = vstary(i-1,j-2)
      u2 = vstarx(i-1,j-1) ; v2 = vstary(i-1,j-1)
      u3 = vstarx(i-1,j-0) ; v3 = vstary(i-1,j-0)
      u4 = vstarx(i-1,j+1) ; v4 = vstary(i-1,j+1)
      call limiter_phi_4(x1,y1,x2,y2,x3,y3,x4,y4,&
            u1,v1,u2,v2,u3,v3,u4,v4,phi3)

      x1 = vertx(i-0,j-2) ; y1 = verty(i-0,j-2)
      x2 = vertx(i-0,j-1) ; y2 = verty(i-0,j-1)
      x3 = vertx(i-0,j-0) ; y3 = verty(i-0,j-0)
      x4 = vertx(i-0,j+1) ; y4 = verty(i-0,j+1)
      u1 = vstarx(i-0,j-2) ; v1 = vstary(i-0,j-2)
      u2 = vstarx(i-0,j-1) ; v2 = vstary(i-0,j-1)
      u3 = vstarx(i-0,j-0) ; v3 = vstary(i-0,j-0)
      u4 = vstarx(i-0,j+1) ; v4 = vstary(i-0,j+1)
      call limiter_phi_4(x1,y1,x2,y2,x3,y3,x4,y4,&
            u1,v1,u2,v2,u3,v3,u4,v4,phi4)

      ! phi(i,j) = (1.d0-phi1) * (1.d0-phi2) * (1.d0-phi3) * (1.d0-phi4) 

      ! phi(i,j) = (1.d0-phi1) + (1.d0-phi2) + (1.d0-phi3) + (1.d0-phi4) 

      phi(i,j) = phi1 + phi2 + phi3 + phi4 

      if( phi(i,j) > 4.d0 - 1.d-2  ) then 
         phi(i,j) = 0.d0 !>  
      else 
         phi(i,j) = 1.d0 
      endif 
   enddo
   enddo 
   return 

   do j = 2 , ny-1
   do i = 1 , 1
      ! x1 = vertx(i-2,j-1) ; y1 = verty(i-2,j-1)
      ! x2 = vertx(i-1,j-1) ; y2 = verty(i-1,j-1)
      ! x3 = vertx(i-0,j-1) ; y3 = verty(i-0,j-1)
      ! x4 = vertx(i+1,j-1) ; y4 = verty(i+1,j-1)
      ! u1 = vertx(i-2,j-1) ; v1 = verty(i-2,j-1)
      ! u2 = vertx(i-1,j-1) ; v2 = verty(i-1,j-1)
      ! u3 = vertx(i-0,j-1) ; v3 = verty(i-0,j-1)
      ! u4 = vertx(i+1,j-1) ; v4 = verty(i+1,j-1)
      ! call limiter_phi_4(x1,y1,x2,y2,x3,y3,x4,y4,&
      !       u1,v1,u2,v2,u3,v3,u4,v4,phi1)

      ! x1 = vertx(i-2,j-0) ; y1 = verty(i-2,j-0)
      ! x2 = vertx(i-1,j-0) ; y2 = verty(i-1,j-0)
      ! x3 = vertx(i-0,j-0) ; y3 = verty(i-0,j-0)
      ! x4 = vertx(i+1,j-0) ; y4 = verty(i+1,j-0)
      ! u1 = vertx(i-2,j-0) ; v1 = verty(i-2,j-0)
      ! u2 = vertx(i-1,j-0) ; v2 = verty(i-1,j-0)
      ! u3 = vertx(i-0,j-0) ; v3 = verty(i-0,j-0)
      ! u4 = vertx(i+1,j-0) ; v4 = verty(i+1,j-0)
      ! call limiter_phi_4(x1,y1,x2,y2,x3,y3,x4,y4,&
      !       u1,v1,u2,v2,u3,v3,u4,v4,phi2)

      x1 = vertx(i-1,j-2) ; y1 = verty(i-1,j-2)
      x2 = vertx(i-1,j-1) ; y2 = verty(i-1,j-1)
      x3 = vertx(i-1,j-0) ; y3 = verty(i-1,j-0)
      x4 = vertx(i-1,j+1) ; y4 = verty(i-1,j+1)
      u1 = vertx(i-1,j-2) ; v1 = verty(i-1,j-2)
      u2 = vertx(i-1,j-1) ; v2 = verty(i-1,j-1)
      u3 = vertx(i-1,j-0) ; v3 = verty(i-1,j-0)
      u4 = vertx(i-1,j+1) ; v4 = verty(i-1,j+1)
      call limiter_phi_4(x1,y1,x2,y2,x3,y3,x4,y4,&
            u1,v1,u2,v2,u3,v3,u4,v4,phi3)

      x1 = vertx(i-0,j-2) ; y1 = verty(i-0,j-2)
      x2 = vertx(i-0,j-1) ; y2 = verty(i-0,j-1)
      x3 = vertx(i-0,j-0) ; y3 = verty(i-0,j-0)
      x4 = vertx(i-0,j+1) ; y4 = verty(i-0,j+1)
      u1 = vertx(i-0,j-2) ; v1 = verty(i-0,j-2)
      u2 = vertx(i-0,j-1) ; v2 = verty(i-0,j-1)
      u3 = vertx(i-0,j-0) ; v3 = verty(i-0,j-0)
      u4 = vertx(i-0,j+1) ; v4 = verty(i-0,j+1)
      call limiter_phi_4(x1,y1,x2,y2,x3,y3,x4,y4,&
            u1,v1,u2,v2,u3,v3,u4,v4,phi4)

      ! phi(i,j) = (1.d0-phi3) * (1.d0-phi4) 

      phi(i,j) = (1.d0-phi3) + (1.d0-phi4) 

      phi(i,j) = 0.d0 
   enddo
   enddo 

   do j = 2 , ny-1
   do i = nx , nx
      ! x1 = vertx(i-2,j-1) ; y1 = verty(i-2,j-1)
      ! x2 = vertx(i-1,j-1) ; y2 = verty(i-1,j-1)
      ! x3 = vertx(i-0,j-1) ; y3 = verty(i-0,j-1)
      ! x4 = vertx(i+1,j-1) ; y4 = verty(i+1,j-1)
      ! u1 = vertx(i-2,j-1) ; v1 = verty(i-2,j-1)
      ! u2 = vertx(i-1,j-1) ; v2 = verty(i-1,j-1)
      ! u3 = vertx(i-0,j-1) ; v3 = verty(i-0,j-1)
      ! u4 = vertx(i+1,j-1) ; v4 = verty(i+1,j-1)
      ! call limiter_phi_4(x1,y1,x2,y2,x3,y3,x4,y4,&
      !       u1,v1,u2,v2,u3,v3,u4,v4,phi1)

      ! x1 = vertx(i-2,j-0) ; y1 = verty(i-2,j-0)
      ! x2 = vertx(i-1,j-0) ; y2 = verty(i-1,j-0)
      ! x3 = vertx(i-0,j-0) ; y3 = verty(i-0,j-0)
      ! x4 = vertx(i+1,j-0) ; y4 = verty(i+1,j-0)
      ! u1 = vertx(i-2,j-0) ; v1 = verty(i-2,j-0)
      ! u2 = vertx(i-1,j-0) ; v2 = verty(i-1,j-0)
      ! u3 = vertx(i-0,j-0) ; v3 = verty(i-0,j-0)
      ! u4 = vertx(i+1,j-0) ; v4 = verty(i+1,j-0)
      ! call limiter_phi_4(x1,y1,x2,y2,x3,y3,x4,y4,&
      !       u1,v1,u2,v2,u3,v3,u4,v4,phi2)

      x1 = vertx(i-1,j-2) ; y1 = verty(i-1,j-2)
      x2 = vertx(i-1,j-1) ; y2 = verty(i-1,j-1)
      x3 = vertx(i-1,j-0) ; y3 = verty(i-1,j-0)
      x4 = vertx(i-1,j+1) ; y4 = verty(i-1,j+1)
      u1 = vertx(i-1,j-2) ; v1 = verty(i-1,j-2)
      u2 = vertx(i-1,j-1) ; v2 = verty(i-1,j-1)
      u3 = vertx(i-1,j-0) ; v3 = verty(i-1,j-0)
      u4 = vertx(i-1,j+1) ; v4 = verty(i-1,j+1)
      call limiter_phi_4(x1,y1,x2,y2,x3,y3,x4,y4,&
            u1,v1,u2,v2,u3,v3,u4,v4,phi3)

      x1 = vertx(i-0,j-2) ; y1 = verty(i-0,j-2)
      x2 = vertx(i-0,j-1) ; y2 = verty(i-0,j-1)
      x3 = vertx(i-0,j-0) ; y3 = verty(i-0,j-0)
      x4 = vertx(i-0,j+1) ; y4 = verty(i-0,j+1)
      u1 = vertx(i-0,j-2) ; v1 = verty(i-0,j-2)
      u2 = vertx(i-0,j-1) ; v2 = verty(i-0,j-1)
      u3 = vertx(i-0,j-0) ; v3 = verty(i-0,j-0)
      u4 = vertx(i-0,j+1) ; v4 = verty(i-0,j+1)
      call limiter_phi_4(x1,y1,x2,y2,x3,y3,x4,y4,&
            u1,v1,u2,v2,u3,v3,u4,v4,phi4)

      ! phi(i,j) = (1.d0-phi3) * (1.d0-phi4) 

      phi(i,j) = (1.d0-phi3) + (1.d0-phi4) 

      phi(i,j) = 0.d0 

   enddo
   enddo 

   do j = 1 , 1
   do i = 2 , nx-1
      x1 = vertx(i-2,j-1) ; y1 = verty(i-2,j-1)
      x2 = vertx(i-1,j-1) ; y2 = verty(i-1,j-1)
      x3 = vertx(i-0,j-1) ; y3 = verty(i-0,j-1)
      x4 = vertx(i+1,j-1) ; y4 = verty(i+1,j-1)
      u1 = vertx(i-2,j-1) ; v1 = verty(i-2,j-1)
      u2 = vertx(i-1,j-1) ; v2 = verty(i-1,j-1)
      u3 = vertx(i-0,j-1) ; v3 = verty(i-0,j-1)
      u4 = vertx(i+1,j-1) ; v4 = verty(i+1,j-1)
      call limiter_phi_4(x1,y1,x2,y2,x3,y3,x4,y4,&
            u1,v1,u2,v2,u3,v3,u4,v4,phi1)

      x1 = vertx(i-2,j-0) ; y1 = verty(i-2,j-0)
      x2 = vertx(i-1,j-0) ; y2 = verty(i-1,j-0)
      x3 = vertx(i-0,j-0) ; y3 = verty(i-0,j-0)
      x4 = vertx(i+1,j-0) ; y4 = verty(i+1,j-0)
      u1 = vertx(i-2,j-0) ; v1 = verty(i-2,j-0)
      u2 = vertx(i-1,j-0) ; v2 = verty(i-1,j-0)
      u3 = vertx(i-0,j-0) ; v3 = verty(i-0,j-0)
      u4 = vertx(i+1,j-0) ; v4 = verty(i+1,j-0)
      call limiter_phi_4(x1,y1,x2,y2,x3,y3,x4,y4,&
            u1,v1,u2,v2,u3,v3,u4,v4,phi2)

      ! x1 = vertx(i-1,j-2) ; y1 = verty(i-1,j-2)
      ! x2 = vertx(i-1,j-1) ; y2 = verty(i-1,j-1)
      ! x3 = vertx(i-1,j-0) ; y3 = verty(i-1,j-0)
      ! x4 = vertx(i-1,j+1) ; y4 = verty(i-1,j+1)
      ! u1 = vertx(i-1,j-2) ; v1 = verty(i-1,j-2)
      ! u2 = vertx(i-1,j-1) ; v2 = verty(i-1,j-1)
      ! u3 = vertx(i-1,j-0) ; v3 = verty(i-1,j-0)
      ! u4 = vertx(i-1,j+1) ; v4 = verty(i-1,j+1)
      ! call limiter_phi_4(x1,y1,x2,y2,x3,y3,x4,y4,&
      !       u1,v1,u2,v2,u3,v3,u4,v4,phi3)

      ! x1 = vertx(i-0,j-2) ; y1 = verty(i-0,j-2)
      ! x2 = vertx(i-0,j-1) ; y2 = verty(i-0,j-1)
      ! x3 = vertx(i-0,j-0) ; y3 = verty(i-0,j-0)
      ! x4 = vertx(i-0,j+1) ; y4 = verty(i-0,j+1)
      ! u1 = vertx(i-0,j-2) ; v1 = verty(i-0,j-2)
      ! u2 = vertx(i-0,j-1) ; v2 = verty(i-0,j-1)
      ! u3 = vertx(i-0,j-0) ; v3 = verty(i-0,j-0)
      ! u4 = vertx(i-0,j+1) ; v4 = verty(i-0,j+1)
      ! call limiter_phi_4(x1,y1,x2,y2,x3,y3,x4,y4,&
      !       u1,v1,u2,v2,u3,v3,u4,v4,phi4)

      ! phi(i,j) = (1.d0-phi1) * (1.d0-phi2) 

      phi(i,j) = (1.d0-phi1) + (1.d0-phi2) 

      phi(i,j) = 0.d0 

   enddo
   enddo 

   do j = ny , ny
   do i = 2 , nx-1
      x1 = vertx(i-2,j-1) ; y1 = verty(i-2,j-1)
      x2 = vertx(i-1,j-1) ; y2 = verty(i-1,j-1)
      x3 = vertx(i-0,j-1) ; y3 = verty(i-0,j-1)
      x4 = vertx(i+1,j-1) ; y4 = verty(i+1,j-1)
      u1 = vertx(i-2,j-1) ; v1 = verty(i-2,j-1)
      u2 = vertx(i-1,j-1) ; v2 = verty(i-1,j-1)
      u3 = vertx(i-0,j-1) ; v3 = verty(i-0,j-1)
      u4 = vertx(i+1,j-1) ; v4 = verty(i+1,j-1)
      call limiter_phi_4(x1,y1,x2,y2,x3,y3,x4,y4,&
            u1,v1,u2,v2,u3,v3,u4,v4,phi1)

      x1 = vertx(i-2,j-0) ; y1 = verty(i-2,j-0)
      x2 = vertx(i-1,j-0) ; y2 = verty(i-1,j-0)
      x3 = vertx(i-0,j-0) ; y3 = verty(i-0,j-0)
      x4 = vertx(i+1,j-0) ; y4 = verty(i+1,j-0)
      u1 = vertx(i-2,j-0) ; v1 = verty(i-2,j-0)
      u2 = vertx(i-1,j-0) ; v2 = verty(i-1,j-0)
      u3 = vertx(i-0,j-0) ; v3 = verty(i-0,j-0)
      u4 = vertx(i+1,j-0) ; v4 = verty(i+1,j-0)
      call limiter_phi_4(x1,y1,x2,y2,x3,y3,x4,y4,&
            u1,v1,u2,v2,u3,v3,u4,v4,phi2)

      ! x1 = vertx(i-1,j-2) ; y1 = verty(i-1,j-2)
      ! x2 = vertx(i-1,j-1) ; y2 = verty(i-1,j-1)
      ! x3 = vertx(i-1,j-0) ; y3 = verty(i-1,j-0)
      ! x4 = vertx(i-1,j+1) ; y4 = verty(i-1,j+1)
      ! u1 = vertx(i-1,j-2) ; v1 = verty(i-1,j-2)
      ! u2 = vertx(i-1,j-1) ; v2 = verty(i-1,j-1)
      ! u3 = vertx(i-1,j-0) ; v3 = verty(i-1,j-0)
      ! u4 = vertx(i-1,j+1) ; v4 = verty(i-1,j+1)
      ! call limiter_phi_4(x1,y1,x2,y2,x3,y3,x4,y4,&
      !       u1,v1,u2,v2,u3,v3,u4,v4,phi3)

      ! x1 = vertx(i-0,j-2) ; y1 = verty(i-0,j-2)
      ! x2 = vertx(i-0,j-1) ; y2 = verty(i-0,j-1)
      ! x3 = vertx(i-0,j-0) ; y3 = verty(i-0,j-0)
      ! x4 = vertx(i-0,j+1) ; y4 = verty(i-0,j+1)
      ! u1 = vertx(i-0,j-2) ; v1 = verty(i-0,j-2)
      ! u2 = vertx(i-0,j-1) ; v2 = verty(i-0,j-1)
      ! u3 = vertx(i-0,j-0) ; v3 = verty(i-0,j-0)
      ! u4 = vertx(i-0,j+1) ; v4 = verty(i-0,j+1)
      ! call limiter_phi_4(x1,y1,x2,y2,x3,y3,x4,y4,&
      !       u1,v1,u2,v2,u3,v3,u4,v4,phi4)

      ! phi(i,j) = (1.d0-phi1) * (1.d0-phi2) 

      phi(i,j) = (1.d0-phi1) + (1.d0-phi2) 

      phi(i,j) = 0.d0 

   enddo
   enddo 

   do j = 1, ny 
   do i = 1, nx 
      if ( phi(i,j) > 0.5d0 ) then 
         phi(i,j) = 1.d0 
      elseif ( phi(i,j) > -0.5d0 ) then 
         phi(i,j) = 0.d0 
      else 
         phi(i,j) = 0.d0 
      endif 
      ! print*, i,j,  phi(i,j) 
   enddo
   enddo 
   ! pause 1293

end subroutine cal_limiter_phi


subroutine limiter_phi_4(x1,y1,x2,y2,x3,y3,x4,y4,&
            u1,v1,u2,v2,u3,v3,u4,v4,phi)
   implicit none
   
	real*8,intent(in)   :: x1,y1,x2,y2,x3,y3,x4,y4
   real*8,intent(in)   :: u1,v1,u2,v2,u3,v3,u4,v4
   real*8,intent(out)  :: phi 

   real*8  :: delu,delv,delx,dely,del,rrk,rlk 

   delu = (u3-u2) / sqrt((u3-u2)**2+(v3-v2)**2)
   delv = (v3-v2) / sqrt((u3-u2)**2+(v3-v2)**2)
   delx = (x3-x2) / sqrt((x3-x2)**2+(y3-y2)**2)
   dely = (y3-y2) / sqrt((x3-x2)**2+(y3-y2)**2)
   del  = sqrt((u3-u2)**2+(v3-v2)**2) / sqrt((x3-x2)**2+(y3-y2)**2)

   rrk = ((u4-u3)*delu+(v4-v3)*delv) / ((x4-x3)*delx+(y4-y3)*dely) / del 
   rlk = ((u2-u1)*delu+(v2-v1)*delv) / ((x2-x1)*delx+(y2-y1)*dely) / del 

   phi = max(0.d0,min(min( min(0.5d0*(rrk+rlk),2.d0*rlk),2.d0*rrk), 1.d0 ) )

end subroutine limiter_phi_4
