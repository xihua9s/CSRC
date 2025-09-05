!> 20230112
!> xihua
!> compute pre1 and pre1_in
!> pre1:    [1],[2],[3],[4]
!> pre1_in: (1),(2),(3),(4)
! -------------------------
! |          |            |
! |          |            |
! |   [4]   (3)    [3]    |
! |          |            |
! |          |            |
! ----(4)----------(2)-----
! |          |            |
! |          |            |
! |   [1]   (1)    [2]    |
! |          |            |
! |          |            |
! -------------------------
subroutine cal_mjsls_details_subcell_edge_pre4_in_ca(nx,ny,eos,vertx,verty, &
   vstarx,vstary,den,pre,sound,gamar,amass,xc,yc,pre1,pre1_in,smass,smass_tri)
	implicit none 
   integer :: nx,ny
	real*8 :: den(1:nx,1:ny),pre(1:nx,1:ny)
	real*8 :: sound(1:nx,1:ny),amass(1:nx,1:ny),smass(1:nx,1:ny,1:4),smass_tri(1:nx,1:ny,1:8)
	real*8 :: vertx(0:nx,0:ny),verty(0:nx,0:ny)
	real*8 :: vstarx(0:nx,0:ny),vstary(0:nx,0:ny)
   integer :: eos(1:nx,1:ny)
	real*8  :: gamar 
	real*8,intent(out)  :: pre1(4,1:nx,1:ny),pre1_in(4,1:nx,1:ny)
	real*8,intent(out)  :: xc(1:nx,1:ny),yc(1:nx,1:ny)
   
   integer  :: i,j 
   real*8   :: phi(1:nx,1:ny)

   !> compute pre1 and pre1_in
	do j = 1 , ny 
	do i = 1 , nx 
      call cal_eos(eos(i,j),gamar)

		! call cal_subcell_edge_p4_in_ca(                                     &
		! 		vertx(i-1,j-1),verty(i-1,j-1),vertx(i,j-1),verty(i,j-1),     &
		! 		vertx(i,j),verty(i,j),vertx(i-1,j),verty(i-1,j),             &
		! 		vstarx(i-1,j-1),vstary(i-1,j-1),vstarx(i,j-1),vstary(i,j-1), &
		! 		vstarx(i,j),vstary(i,j),vstarx(i-1,j),vstary(i-1,j),         &
      !       xc(i,j),yc(i,j),den(i,j),pre(i,j),sound(i,j),pre1(1:4,i,j),  &
      !       gamar,amass(i,j),pre1_in(1:4,i,j),phi(i,j),smass(i,j,1:4),   &
      !       smass_tri(i,j,1:8) )     
     
	enddo
   enddo

end subroutine cal_mjsls_details_subcell_edge_pre4_in_ca

 !> 20250724
 !> edge viscosity 
 !> Xihua
subroutine cal_subcell_edge_p4_in_ca(x1,y1,x2,y2,x3,y3,x4,y4,&
   u1,v1,u2,v2,u3,v3,u4,v4,xc,yc, &
   den,pre,sound,pre1,gamar,amass,pre1_in,phi)
	implicit none 
	real*8,intent(in)   :: x1,y1,x2,y2,x3,y3,x4,y4
	real*8,intent(in)   :: u1,v1,u2,v2,u3,v3,u4,v4
	real*8,intent(in)   :: den,pre,sound,gamar,amass,phi
	real*8,intent(out)  :: pre1(4),pre1_in(4),xc,yc

	integer  :: i,i1,i2,side,n
	real*8   :: x12,y12,x23,y23,x34,y34,x41,y41,sub_area(4),cy_volm
	real*8   :: xa(4),ya(4),ua(4),va(4),alength,anormalx,anormaly,ca_area
	real*8   :: delta_v,clength,area,xx1,yy1,hbar,delta_tau,xx2,yy2
   real*8   :: ustar,vstar,hbar4(4),pstar4(4),rwai,rnei,uc,vc,tden 
   real*8   :: A12(2,4),A2c(2,4),Ac1(2,4),D12(4),AAA,BBB,kkk
   real*8   :: uu1,vv1,uu2,vv2,BB(4)

	xc = (x1+x2+x3+x4) / 4.d0 
	yc = (y1+y2+y3+y4) / 4.d0 

   ! call cal_mass_center(x1,y1,x2,y2,x3,y3,x4,y4,xc,yc)

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
      area = area + ca_area(xx1,yy1,xx2,yy2)      
      clength = clength + alength
   enddo 
   hbar = area/clength * 2.d0 * 2.d0 
   delta_tau = delta_v * hbar / area / den / sound

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
      uu1 = ua(i1) ; vv1 = va(i1)
      uu2 = ua(i2) ; vv2 = va(i2)
      call triangle_geo_f1(xx1,yy1,xx2,yy2,xc,yc,A12(1:2,i),A2c(1:2,i),Ac1(1:2,i))
      D12(i) = uu1*A2c(1,i)+vv1*A2c(2,i) + uu2*Ac1(1,i)+vv2*Ac1(2,i)
      call triangle_geo_f2(xx1,yy1,xx2,yy2,xc,yc,uu1,vv1,uu2,vv2,BB(i))
   enddo 

   AAA = den*sound*( uc*(BB(1)*A12(1,1)+BB(2)*A12(1,2)+BB(3)*A12(1,3)+BB(4)*A12(1,4)) + &
                     vc*(BB(1)*A12(2,1)+BB(2)*A12(2,2)+BB(3)*A12(2,3)+BB(4)*A12(2,4)) )
   BBB = pre*sum(BB(:)) - den*sound*dot_product(D12(:),BB(:))

   if ( abs(AAA) < 1.d-16 ) then 
      kkk = 0.d0 
      do i = 1 , 4
         pre1(i) = pre !-den*sound*kkk*(uc*A12(1,i)+vc*A12(2,i)) + pre - den*sound*D12(i)
      enddo       
   else 
      kkk = BBB / AAA
      do i = 1 , 4
         pre1(i) = -den*sound*kkk*(uc*A12(1,i)+vc*A12(2,i)) + pre - den*sound*D12(i) &
                  + (gamar+1.d0)/2.d0*den**3*sound**2*delta_tau**2
      enddo 
   endif 

    if( delta_tau < 0.d0 ) then 
      ! pre1 = pre - (den*sound)**2*delta_tau  &
      !    + (gamar+1.d0)/2.d0*den**3*sound**2*delta_tau**2    
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
 
   !> second method
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
			area = area + ca_area(xx1,yy1,xx2,yy2)         
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
      pre1_in(side)  = sound**2 * (amass/2.d0/(sub_area(i1)+sub_area(i2)) - den) !* 0.d0 
   enddo 

end subroutine cal_subcell_edge_p4_in_ca


