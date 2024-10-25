

subroutine Update_Mesh_Phys_Cat_Xihua_0_1017(den,velx,vely,pre,sound,tau,enin,Ener,EOS, &
            Vstarx,Vstary,Pstar1,mass,vertx,verty,area,volm,nx,ny,dt)
   use COMMON_TYPE_EOS
   Use COMMON_XXX_Constant
   implicit none

   real*8,intent(inout) :: den(1:nx,1:ny),velx(1:nx,1:ny),vely(1:nx,1:ny)
   real*8,intent(inout) :: pre(1:nx,1:ny),mass(1:nx,1:ny)
   real*8,intent(inout) :: sound(1:nx,1:ny),tau(1:nx,1:ny)
   real*8,intent(inout) :: enin(1:nx,1:ny),Ener(1:nx,1:ny)
   integer(KNDI),intent(in) :: EOS(1:nx,1:ny)
   real*8,intent(inout) :: vertx(0:nx,0:ny),verty(0:nx,0:ny)
   real*8,intent(inout) :: area(1:nx,1:ny),volm(1:nx,1:ny)
   real*8,intent(in)  :: Vstarx(0:nx,0:ny),Vstary(0:nx,0:ny)
   real*8,intent(in)  :: Pstar1(0:nx,0:ny)
   integer,intent(in) :: nx,ny
   real*8,intent(in)  :: dt
   
   real*8 :: tempttau(1:nx,1:ny),temptvelx(1:nx,1:ny)
   real*8 :: temptEner(1:nx,1:ny),temptvely(1:nx,1:ny),temptein(1:nx,1:ny)
    
   integer  :: i,j,ir,jr,ir1,jr1,ix,iy,side,nn
   real*8   :: sum1,sum2x,sum2y,sum3,xr2,yr2,cy_area,cy_volm,sum4,epseps,v1,v2,ratiox,ratioy
   real*8   :: xr,yr,xr1,yr1,rr1,rr2,alengthn,anormalxn,anormalyn,vv,tmp_velx,tmp_vely,Sum5
   real*8   :: vbarx(4),vbary(4),vbarx1(4),vbary1(4),v0
   real*8   :: avel(nx),rr,xmid1,ymid1,xmid2,ymid2,rh,ttenin,sigma
   
   nn = 10000
   sigma = 1.d-10
   do i=2,nx
   do j=1,ny

      Sum1=0.d0
      Sum2x=0.d0
      Sum2y=0.d0
      Sum3=0.d0
      Sum4=0.d0
      Sum5=0.d0
      
      tmp_velx = 0.d0
      tmp_vely = 0.d0
      
      rr = (vertx(i,j-1)-vertx(i-1,j-1))**2 + (verty(i,j-1)-verty(i-1,j-1))**2 
      xmid1 = (vertx(i-1,j)+vertx(i-1,j-1))/2.d0
      ymid1 = (verty(i-1,j)+verty(i-1,j-1))/2.d0
      xmid2 = (vertx(i,j)+vertx(i,j-1))/2.d0
      ymid2 = (verty(i,j)+verty(i,j-1))/2.d0
      rh = ( (xmid2-xmid1)**2 + (ymid2-ymid1)**2 )
      
      do side = 1 , 4    !> 4 means quadrangle
         select case(side)
         case(1)
            ir = i-1
            jr = j-1
            ir1 = i  
            jr1 = j-1
         case(2)
            ir = i
            jr = j-1
            ir1 = i
            jr1 = j
         case(3)
            ir = i
            jr = j
            ir1 = i-1
            jr1 = j  
         case(4)
            ir=i-1
            jr=j
            ir1=i-1
            jr1=j-1         
         case default
            stop 'side is Error, pls Call Xihua'
         end select
         
         tmp_velx = tmp_velx + sigma * Vstarx(ir,jr) + sigma * Vstarx(ir1,jr1)
         tmp_vely = tmp_vely + sigma * Vstary(ir,jr) + sigma * Vstary(ir1,jr1)         
         
      enddo
      
      tmp_velx = tmp_velx / 4.d0 + (1.d0-2.d0*sigma)*velx(i,j)*sqrt(rh/rr) 
      tmp_vely = tmp_vely / 4.d0 + (1.d0-2.d0*sigma)*vely(i,j)*sqrt(rh/rr) 
      
      do side = 1 , 4    !> 4 means quadrangle
         select case(side)
         case(1)
            ir = i-1
            jr = j-1
            ir1 = i  
            jr1 = j-1
         case(2)
            ir = i
            jr = j-1
            ir1 = i
            jr1 = j
         case(3)
            ir = i
            jr = j
            ir1 = i-1
            jr1 = j  
         case(4)
            ir=i-1
            jr=j
            ir1=i-1
            jr1=j-1         
         case default
            stop 'side is Error, pls Call Xihua'
         end select

         xr=vertx(ir,jr)           
         yr=verty(ir,jr)
         xr1=vertx(ir1,jr1)
         yr1=verty(ir1,jr1)
      
         call Length(xr,yr,xr1,yr1,alengthn)
         alengthn = alengthn
         call Normal(xr,yr,xr1,yr1,anormalxn,anormalyn)
         anormalxn=-anormalxn
         anormalyn=-anormalyn

         Sum1 = Sum1 + anormalxn*alengthn*(Vstarx(ir,jr)+Vstarx(ir1,jr1))/2.d0 + &
                       anormalyn*alengthn*(Vstary(ir,jr)+Vstary(ir1,jr1))/2.d0

         Sum2x=Sum2x+anormalxn*alengthn*(Pstar1(ir,jr)+Pstar1(ir1,jr1))/2.d0
         Sum2y=Sum2y+anormalyn*alengthn*(Pstar1(ir,jr)+Pstar1(ir1,jr1))/2.d0
         
         Sum3 = Sum3 + anormalxn*alengthn*(Pstar1(ir,jr)*Vstarx(ir,jr) + &
                                  Pstar1(ir1,jr1)*Vstarx(ir1,jr1))/2.d0 + &
                       anormalyn*alengthn*(Pstar1(ir,jr)*Vstary(ir,jr) + &
                                  Pstar1(ir1,jr1)*Vstary(ir1,jr1))/2.d0         

         Sum4 = Sum4 + anormalxn*alengthn*(Pstar1(ir,jr)*(Vstarx(ir,jr)-tmp_velx) + &
                      Pstar1(ir1,jr1) * (Vstarx(ir1,jr1)- tmp_velx)     )/2.d0 + &
                       anormalyn*alengthn*(Pstar1(ir,jr)*(Vstary(ir,jr)-tmp_vely) + &
                      Pstar1(ir1,jr1) * (Vstary(ir1,jr1)- tmp_vely)    )/2.d0

      enddo

      tempttau(i,j)  = tau(i,j)  + Sum1  * dt/mass(i,j)
      temptvelx(i,j) = velx(i,j) - Sum2x * dt/mass(i,j)
      temptvely(i,j) = vely(i,j) - Sum2y * dt/mass(i,j)
      temptEner(i,j) = Ener(i,j) - Sum3  * dt/mass(i,j)      
      temptein(i,j) = enin(i,j)  - Sum4  * dt/mass(i,j)
      
   enddo 
   enddo


   do i=1,1
   do j=1,ny

      Sum1=0.d0
      Sum2x=0.d0
      Sum2y=0.d0
      Sum3=0.d0
      Sum4=0.d0
      Sum5=0.d0
      
      tmp_velx = 0.d0
      tmp_vely = 0.d0
      
      rr = (vertx(i,j-1)-vertx(i-1,j-1))**2 + (verty(i,j-1)-verty(i-1,j-1))**2 
      xmid1 = (vertx(i-1,j)+vertx(i-1,j-1))/2.d0
      ymid1 = (verty(i-1,j)+verty(i-1,j-1))/2.d0
      xmid2 = (vertx(i,j)+vertx(i,j-1))/2.d0
      ymid2 = (verty(i,j)+verty(i,j-1))/2.d0
      rh = ( (xmid2-xmid1)**2 + (ymid2-ymid1)**2 )
      
      do side = 1 , 3    !> 3 means triangle
         select case(side)
         case(1)
            ir = i-1
            jr = j-1
            ir1 = i  
            jr1 = j-1
         case(2)
            ir = i
            jr = j-1
            ir1 = i
            jr1 = j
         case(3)
            ir = i
            jr = j
            ir1 = i-1
            jr1 = j  
         ! case(4)
            ! ir=i-1
            ! jr=j
            ! ir1=i-1
            ! jr1=j-1         
         case default
            stop 'side is Error, pls Call Xihua'
         end select
                           
         tmp_velx = tmp_velx + sigma * Vstarx(ir,jr) + sigma * Vstarx(ir1,jr1)
         tmp_vely = tmp_vely + sigma * Vstary(ir,jr) + sigma * Vstary(ir1,jr1)         
         
      enddo
      
      tmp_velx = tmp_velx / 3.d0 + (1.d0-2.d0*sigma)*velx(i,j)*sqrt(rh/rr) 
      tmp_vely = tmp_vely / 3.d0 + (1.d0-2.d0*sigma)*vely(i,j)*sqrt(rh/rr) 
      
      do side = 1 , 3    !> 3 means triangle
         select case(side)
         case(1)
            ir = i-1
            jr = j-1
            ir1 = i  
            jr1 = j-1
         case(2)
            ir = i
            jr = j-1
            ir1 = i
            jr1 = j
         case(3)
            ir = i
            jr = j
            ir1 = i-1
            jr1 = j  
         ! case(4)
            ! ir=i-1
            ! jr=j
            ! ir1=i-1
            ! jr1=j-1         
         case default
            stop 'side is Error, pls Call Xihua'
         end select

         xr=vertx(ir,jr)           
         yr=verty(ir,jr)
         xr1=vertx(ir1,jr1)
         yr1=verty(ir1,jr1)
      
         call Length(xr,yr,xr1,yr1,alengthn)
         alengthn = alengthn
         call Normal(xr,yr,xr1,yr1,anormalxn,anormalyn)
         anormalxn=-anormalxn
         anormalyn=-anormalyn
         
         Sum1 = Sum1 + ( anormalxn*alengthn*Vstarx(ir,jr)   +&
                         anormalyn*alengthn*Vstary(ir,jr)       )/2.d0 + &
                       ( anormalxn*alengthn*Vstarx(ir1,jr1) + &
                         anormalyn*alengthn*Vstary(ir1,jr1)     )/2.d0
                         
         Sum2x=Sum2x+anormalxn*alengthn*(Pstar1(ir,jr)+Pstar1(ir1,jr1))/2.d0
         Sum2y=Sum2y+anormalyn*alengthn*(Pstar1(ir,jr)+Pstar1(ir1,jr1))/2.d0
         
         Sum3 = Sum3 + anormalxn*alengthn*(Pstar1(ir,jr)*Vstarx(ir,jr) + &
                                            Pstar1(ir1,jr1)*Vstarx(ir1,jr1))/2.d0 + &
                        anormalyn*alengthn*(Pstar1(ir,jr)*Vstary(ir,jr) + &
                                            Pstar1(ir1,jr1)*Vstary(ir1,jr1))/2.d0

         Sum4 = Sum4 + ( anormalxn*alengthn*Pstar1(ir,jr)*Vstarx(ir,jr)     - &
                         anormalxn*alengthn*Pstar1(ir,jr)*tmp_velx          + &
                         anormalxn*alengthn*Pstar1(ir1,jr1)*Vstarx(ir1,jr1) - &
                         anormalxn*alengthn*Pstar1(ir1,jr1)*tmp_velx            )/2.d0 + &
                       ( anormalyn*alengthn*Pstar1(ir,jr)*Vstary(ir,jr)     - &
                         anormalyn*alengthn*Pstar1(ir,jr)*tmp_vely          + &
                         anormalyn*alengthn*Pstar1(ir1,jr1)*Vstary(ir1,jr1) - &
                         anormalyn*alengthn*Pstar1(ir1,jr1)*tmp_vely            )/2.d0 
         
      enddo

      tempttau(i,j)  = tau(i,j)  + Sum1  * dt/mass(i,j)
      temptvelx(i,j) = velx(i,j) - Sum2x * dt/mass(i,j)
      temptvely(i,j) = vely(i,j) - Sum2y * dt/mass(i,j)
      temptEner(i,j) = Ener(i,j) - Sum3  * dt/mass(i,j)      
      temptein(i,j) = enin(i,j)  - Sum4  * dt/mass(i,j)
   enddo 
   enddo
   
   ttenin = sum(temptein(1,1:ny))/ny
   temptein(1,1:ny) = ttenin

   do i=1,nx
   do j=1,ny  
      
      Ener(i,j) = temptEner(i,j) 
      enin(i,j) = temptein(i,j)
      vv = (Ener(i,j) - enin(i,j))*2.d0

      if(vv < 0.d0) then
         temptvelx(i,j) = -velx(i,j) 
         temptvely(i,j) = -vely(i,j) 
         vv = velx(i,j)*velx(i,j) + vely(i,j)*vely(i,j)
      endif      
      
      tau(i,j)  = tempttau(i,j)
      velx(i,j) = temptvelx(i,j)
      vely(i,j) = temptvely(i,j)
      den(i,j)  = 1.d0/tau(i,j)
      !enin(i,j) = Ener(i,j) - (velx(i,j)**2 + vely(i,j)**2)/2.d0      

      v2 = velx(i,j)*velx(i,j) + vely(i,j)*vely(i,j)
      v1 = temptvelx(i,j)*temptvelx(i,j) + temptvely(i,j)*temptvely(i,j)      
      if(vv<eps) then    
         velx(i,j) = 0.d0
         vely(i,j) = 0.d0
      else
         velx(i,j) = temptvelx(i,j) / dsqrt(v1/vv) 
         vely(i,j) = temptvely(i,j) / dsqrt(v1/vv) 
      endif
      
      Call General_EOS(EOS(i,j))
      pre(i,j)  = enin(i,j)*(gammar-1.d0)*den(i,j)
      sound(i,j)= sqrt(gammar*pre(i,j)/den(i,j))

   enddo 
   enddo 
   
   do i=0,nx
   do j=0,ny
      vertx(i,j) = vertx(i,j) + dt*Vstarx(i,j)
      verty(i,j) = verty(i,j) + dt*Vstary(i,j)
   enddo
   enddo 

end subroutine Update_Mesh_Phys_Cat_Xihua_0_1017

