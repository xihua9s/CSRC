

subroutine update_vertex(nx,ny,vertx,verty,vertx_new,verty_new,velx,vely,dt)
   implicit none 
   integer  :: nx,ny 
	real*8 :: vertx(0:nx,0:ny),verty(0:nx,0:ny),dt 
	real*8 :: vertx_new(0:nx,0:ny),verty_new(0:nx,0:ny)
	real*8 :: velx(0:nx,0:ny),vely(0:nx,0:ny)
   
   integer :: i,j 
   
	do j = 0 , ny 
	do i = 0 , nx 
		vertx_new(i,j) = vertx(i,j) + velx(i,j)*dt
		verty_new(i,j) = verty(i,j) + vely(i,j)*dt
	enddo 
	enddo 

end subroutine update_vertex


subroutine update_den_enin(nx,ny,eos,vertx,verty,vertx_new,verty_new,velx,vely,dt,&
            pre_s,den,enin,amass,gamar,den_new,enin_new,pre_new,sound_new,g)
   implicit none 
   integer :: nx,ny
	real*8 :: den(1:nx,1:ny),amass(1:nx,1:ny),enin(1:nx,1:ny)
	real*8 :: vertx(0:nx,0:ny),verty(0:nx,0:ny)
	real*8 :: vertx_new(0:nx,0:ny),verty_new(0:nx,0:ny)
	real*8 :: velx(0:nx,0:ny),vely(0:nx,0:ny),pre_s(4,1:nx,1:ny)
	real*8  :: gamar,dt,g 
	real*8 :: den_new(1:nx,1:ny),pre_new(1:nx,1:ny)
   real*8 :: sound_new(1:nx,1:ny),enin_new(1:nx,1:ny)
	integer :: eos(1:nx,1:ny)

	integer :: i,j,side,i1,j1,i2,j2
   real*8 :: l4(4),normx(4),normy(4),uv(4),area,x1,x2,y1,y2
   
	do j = 1 , ny
	do i = 1 , nx 
      area=0.d0
		do side = 1 , 4
			select case (side)
			case(1)
				i1 = i-1 ; j1 = j-1
				i2 = i   ; j2 = j-1
			case(2)
				i1 = i   ; j1 = j-1
				i2 = i   ; j2 = j
			case(3)
				i1 = i   ; j1 = j
				i2 = i-1 ; j2 = j
			case(4)
				i1 = i-1 ; j1 = j
				i2 = i-1 ; j2 = j-1
			case default
				stop 'Error,fia'
			end select
			x1=vertx(i1,j1);y1=verty(i1,j1)
			x2=vertx(i2,j2);y2=verty(i2,j2)	            
			call cal_length_normal(x1,x2,y1,y2,l4(side),normx(side),normy(side))
			x1=vertx_new(i1,j1);y1=verty_new(i1,j1)
			x2=vertx_new(i2,j2);y2=verty_new(i2,j2)	            
         area=area+0.5d0*(x1*y2-x2*y1)
		enddo 
		uv(1) = 0.5d0*(l4(4)*normx(4)+l4(1)*normx(1))*velx(i-1,j-1) + &
			     0.5d0*(l4(4)*normy(4)+l4(1)*normy(1))*vely(i-1,j-1)
		uv(2) = 0.5d0*(l4(1)*normx(1)+l4(2)*normx(2))*velx(i-0,j-1) + &
			     0.5d0*(l4(1)*normy(1)+l4(2)*normy(2))*vely(i-0,j-1)
		uv(3) = 0.5d0*(l4(2)*normx(2)+l4(3)*normx(3))*velx(i-0,j-0) + &
			     0.5d0*(l4(2)*normy(2)+l4(3)*normy(3))*vely(i-0,j-0)  
		uv(4) = 0.5d0*(l4(3)*normx(3)+l4(4)*normx(4))*velx(i-1,j-0) + &
			     0.5d0*(l4(3)*normy(3)+l4(4)*normy(4))*vely(i-1,j-0)  
		call cal_eos(eos(i,j),gamar)
		den_new(i,j)   = amass(i,j) /area
		enin_new(i,j)  = enin(i,j) - dot_product(pre_s(1:4,i,j),uv(1:4))*dt / amass(i,j) &
		+ amass(i,j)/4.d0*g*(vely(i-1,j-1)+vely(i-0,j-1)+vely(i-0,j-0)+vely(i-1,j-0))*dt/amass(i,j)
		pre_new(i,j)   = (gamar-1.d0)*den_new(i,j)*enin_new(i,j)
		sound_new(i,j) = sqrt(gamar*pre_new(i,j)/den_new(i,j))

	enddo
	enddo 


end subroutine update_den_enin
