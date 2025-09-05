!> Last Modified : 二  9/26 21:44:20 2023
!> 20230112
!> xihua
!> update the vertex 

module mdu_update
	use fluid_euler_eqn_basic

	implicit none 

contains

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


!> 20230112
!> xihua
!> update the density and internal energy 

subroutine update_den_enin_ca(nx,ny,eos,vertx,verty,vertx_new,verty_new,velx,vely,dt,&
            pre_s,den,enin,amass,gamar,den_new,enin_new,pre_new,sound_new,g)
   implicit none 
   integer :: nx,ny
	real*8  :: den(1:nx,1:ny),amass(1:nx,1:ny),enin(1:nx,1:ny)
	real*8  :: vertx(0:nx,0:ny),verty(0:nx,0:ny)
	real*8  :: vertx_new(0:nx,0:ny),verty_new(0:nx,0:ny)
	real*8  :: velx(0:nx,0:ny),vely(0:nx,0:ny),pre_s(4,1:nx,1:ny)
	real*8  :: gamar,dt,g,cy_volm
	real*8  :: den_new(1:nx,1:ny),pre_new(1:nx,1:ny)
   real*8  :: sound_new(1:nx,1:ny),enin_new(1:nx,1:ny)
	integer :: eos(1:nx,1:ny)

	integer :: i,j,side,i1,j1,i2,j2
   real*8  :: l4(4),normx(4),normy(4),uv(4),area,x1,x2,y1,y2
	real*8  :: subcell_area_ca4(4),subcell_volm_cy_aw4(4)
	real*8  :: subcell_area_cy_aw4(4),ca_area,pre_sr(4),subcell_area_tri(8)
   
	do i = 1 , nx 
	do j = 1 , ny
		area = 0.d0 
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
         area=area+ca_area(x1,y1,x2,y2)
		enddo 

      !> 计算子网格的平几何体积和等效质量
      call cal_subcell_area_ca(vertx_new(i-1,j-1),verty_new(i-1,j-1), &
            vertx_new(i,j-1),verty_new(i,j-1),vertx_new(i,j),verty_new(i,j),  & 
            vertx_new(i-1,j),verty_new(i-1,j),subcell_area_ca4(:),subcell_area_tri(:))  

		uv(1) =( normx(1)*(velx(i-1,j-1)+velx(i-0,j-1))/2.d0 + &
					normy(1)*(vely(i-1,j-1)+vely(i-0,j-1))/2.d0 ) * l4(1)
		uv(2) =( normx(2)*(velx(i-0,j-1)+velx(i-0,j-0))/2.d0 + &
					normy(2)*(vely(i-0,j-1)+vely(i-0,j-0))/2.d0 ) * l4(2)
		uv(3) =( normx(3)*(velx(i-0,j-0)+velx(i-1,j-0))/2.d0 + &
					normy(3)*(vely(i-0,j-0)+vely(i-1,j-0))/2.d0 ) * l4(3)
		uv(4) =( normx(4)*(velx(i-1,j-0)+velx(i-1,j-1))/2.d0 + &
					normy(4)*(vely(i-1,j-0)+vely(i-1,j-1))/2.d0 ) * l4(4)

		call cal_eos(eos(i,j),gamar)
		den_new(i,j)   = amass(i,j) / sum(subcell_area_ca4(:))
		enin_new(i,j)  = enin(i,j) - dot_product(pre_s(1:4,i,j),uv(1:4))*dt / amass(i,j)
		!pre_new(i,j)   = (gamar-1.d0)*den_new(i,j)*enin_new(i,j)
		!sound_new(i,j) = sqrt(gamar*pre_new(i,j)/den_new(i,j))
		call den_ein_2_pre_sos(eos(i,j),den_new(i,j),enin_new(i,j),pre_new(i,j),sound_new(i,j))

	enddo
	enddo 
	
end subroutine update_den_enin_ca


!> 20230112
!> xihua
!> update the density and internal energy 

subroutine update_den_enin_cy(nx,ny,eos,vertx,verty,vertx_new,verty_new,velx,vely,dt,&
            pre_s,den,enin,amass,gamar,den_new,enin_new,pre_new,sound_new,g)
   implicit none 
   integer :: nx,ny
	real*8  :: den(1:nx,1:ny),amass(1:nx,1:ny),enin(1:nx,1:ny)
	real*8  :: vertx(0:nx,0:ny),verty(0:nx,0:ny)
	real*8  :: vertx_new(0:nx,0:ny),verty_new(0:nx,0:ny)
	real*8  :: velx(0:nx,0:ny),vely(0:nx,0:ny),pre_s(4,1:nx,1:ny)
	real*8  :: gamar,dt,g,cy_volm
	real*8  :: den_new(1:nx,1:ny),pre_new(1:nx,1:ny)
   real*8  :: sound_new(1:nx,1:ny),enin_new(1:nx,1:ny)
	integer :: eos(1:nx,1:ny)

	integer :: i,j,side,i1,j1,i2,j2
   real*8  :: l4(4),normx(4),normy(4),uv(4),area,x1,x2,y1,y2
	real*8  :: subcell_area_ca4(4),subcell_volm_cy_aw4(4)
	real*8  :: subcell_area_cy_aw4(4),ca_area,pre_sr(4),subcell_area_tri(8)
   
	do i = 1 , nx 
	do j = 1 , ny
		area = 0.d0 
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
         area=area+ca_area(x1,y1,x2,y2)
		enddo 

      !> 计算子网格的平几何体积和等效质量
      call cal_subcell_area_ca(vertx_new(i-1,j-1),verty_new(i-1,j-1), &
            vertx_new(i,j-1),verty_new(i,j-1),vertx_new(i,j),verty_new(i,j),  & 
            vertx_new(i-1,j),verty_new(i-1,j),subcell_area_ca4(:),subcell_area_tri)  
      subcell_area_cy_aw4(:) = 2.d0/3.d0*subcell_area_ca4(:) + &
                               1.d0/12.d0*area
      subcell_volm_cy_aw4(1) = subcell_area_cy_aw4(1) * verty_new(i-1,j-1)
      subcell_volm_cy_aw4(2) = subcell_area_cy_aw4(2) * verty_new(i  ,j-1)
      subcell_volm_cy_aw4(3) = subcell_area_cy_aw4(3) * verty_new(i  ,j  )
      subcell_volm_cy_aw4(4) = subcell_area_cy_aw4(4) * verty_new(i-1,j  )

		uv(1) = 0.5d0*(l4(4)*normx(4)+l4(1)*normx(1))*velx(i-1,j-1) + &
			     0.5d0*(l4(4)*normy(4)+l4(1)*normy(1))*vely(i-1,j-1)
		uv(2) = 0.5d0*(l4(1)*normx(1)+l4(2)*normx(2))*velx(i-0,j-1) + &
			     0.5d0*(l4(1)*normy(1)+l4(2)*normy(2))*vely(i-0,j-1)
		uv(3) = 0.5d0*(l4(2)*normx(2)+l4(3)*normx(3))*velx(i-0,j-0) + &
			     0.5d0*(l4(2)*normy(2)+l4(3)*normy(3))*vely(i-0,j-0)  
		uv(4) = 0.5d0*(l4(3)*normx(3)+l4(4)*normx(4))*velx(i-1,j-0) + &
			     0.5d0*(l4(3)*normy(3)+l4(4)*normy(4))*vely(i-1,j-0)  
				  
		call cal_eos(eos(i,j),gamar)
		den_new(i,j)   = amass(i,j) / sum(subcell_volm_cy_aw4(:))
		pre_sr(1) = pre_s(1,i,j) * verty(i-1,j-1)
		pre_sr(2) = pre_s(2,i,j) * verty(i ,j-1)
		pre_sr(3) = pre_s(3,i,j) * verty(i ,j)
		pre_sr(4) = pre_s(4,i,j) * verty(i-1,j)
		enin_new(i,j)  = enin(i,j) - dot_product(pre_sr(1:4),uv(1:4))*dt / amass(i,j)
		!pre_new(i,j)   = (gamar-1.d0)*den_new(i,j)*enin_new(i,j)
		!sound_new(i,j) = sqrt(gamar*pre_new(i,j)/den_new(i,j))
		call den_ein_2_pre_sos(eos(i,j),den_new(i,j),enin_new(i,j),pre_new(i,j),sound_new(i,j))

	enddo
	enddo 
	
end subroutine update_den_enin_cy


end module mdu_update