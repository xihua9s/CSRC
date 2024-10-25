
Subroutine Get_Dt(Ini_CDT,MeshCY2D,PhysCY2D)
   use COMMON_TYPE_Condition
   use COMMON_TYPE_Mesh
   use COMMON_TYPE_Phys
   use COMMON_XXX_Constant
   implicit none
   
   Type(COMMON_TYPE_CDT),intent(inout) :: Ini_CDT
   Type(Mesh_Cyl_2D),intent(inout) :: MeshCY2D
   Type(Phys_Cyl_2D),intent(inout) :: PhysCY2D

   Call getdt(PhysCY2D%sound,PhysCY2D%mass,PhysCY2D%tau,MeshCY2D%vertex_x, &
         MeshCY2D%vertex_y,PhysCY2D%Vstar_x,PhysCY2D%Vstar_y,Ini_CDT%nx,    &
         Ini_CDT%ny,Ini_CDT%cfl,Ini_CDT%dt)
   
   if( Ini_CDT%time + Ini_CDT%dt > Ini_CDT%Tstop) then
      Ini_CDT%dt = Ini_CDT%Tstop - Ini_CDT%time
   endif
   
   Ini_CDT%time = Ini_CDT%time + Ini_CDT%dt
   Ini_CDT%it = Ini_CDT%it + 1
   write(*,*) Ini_CDT%it,Ini_CDT%dt,Ini_CDT%time

End Subroutine Get_Dt


 subroutine getdt(sound,amass,tau,vertx,verty,Vstarx,Vstary,nx,ny,cfl,dt)
   implicit none 
   
   integer,intent(in) :: nx,ny
   real*8,intent(in) :: cfl
   real*8,intent(in) :: amass(1:nx,1:ny)
   real*8,intent(in) :: sound(1:nx,1:ny),tau(1:nx,1:ny)
   real*8,intent(in) :: vertx(0:nx,0:ny),verty(0:nx,0:ny)
   real*8,intent(in) :: Vstarx(0:nx,0:ny),Vstary(0:nx,0:ny)
   real*8,intent(out) :: dt

   integer :: i,j,ir,jr,ir1,jr1
   real*8 :: CE,CV,CM,tE,tV,alanmda,xm,ym,xn,yn,alength
   real*8 :: dadt,xr,yr,xr1,yr1,ur,vr,ur1,vr1
   
   CE=0.3d0
   CV=0.1d0
   CM=1.01d0
   tE=9999.d0
   tV=9999.d0
   alanmda=9999.d0

   do i=1,nx
   do j=1,ny         
       alanmda=999999.

       xm=vertx(i-1,j-1)
       ym=verty(i-1,j-1)
       xn=vertx(i,j-1)
       yn=verty(i,j-1)
       call Length(xm,ym,xn,yn,alength)
       if(alength > 1.d-10) alanmda=min(alanmda,alength) 

       xm=vertx(i-1,j-1)
       ym=verty(i-1,j-1)
       xn=vertx(i,j)
       yn=verty(i,j)
       call Length(xm,ym,xn,yn,alength)
       if(alength > 1.d-10) alanmda=min(alanmda,alength) 

       xm=vertx(i-1,j-1)
       ym=verty(i-1,j-1)
       xn=vertx(i-1,j)
       yn=verty(i-1,j)
       call Length(xm,ym,xn,yn,alength)
       if(alength > 1.d-10) alanmda=min(alanmda,alength)

       xm=vertx(i,j-1)
       ym=verty(i,j-1)
       xn=vertx(i,j)
       yn=verty(i,j)
       call Length(xm,ym,xn,yn,alength)
       if(alength > 1.d-10) alanmda=min(alanmda,alength)  

       xm=vertx(i,j-1)
       ym=verty(i,j-1)
       xn=vertx(i-1,j)
       yn=verty(i-1,j)
       call Length(xm,ym,xn,yn,alength)
       if(alength > 1.d-10) alanmda=min(alanmda,alength)    
       
       xm=vertx(i,j)
       ym=verty(i,j)
       xn=vertx(i-1,j)
       yn=verty(i-1,j)
       call Length(xm,ym,xn,yn,alength)
       if(alength > 1.d-10) alanmda=min(alanmda,alength)       	

       tE=min(CE*alanmda/sound(i,j),tE)
   enddo            
   enddo 

   do i=1,nx
   do j=1,ny

      dAdt=0.
      !c1	          
      ir=i-1
      jr=j-1
      ir1=i
      jr1=j-1

      xr=vertx(ir,jr)           
      yr=verty(ir,jr)
      xr1=vertx(ir1,jr1)
      yr1=verty(ir1,jr1)

      ur=Vstarx(ir,jr)    ! u star r
      vr=Vstary(ir,jr)    ! v star r
      ur1=Vstarx(ir1,jr1)
      vr1=Vstary(ir1,jr1)

      dAdt=dAdt+(ur*yr1+vr1*xr-ur1*yr-vr*xr1)*0.5
      !c2
      ir=i
      jr=j-1
      ir1=i
      jr1=j

      xr=vertx(ir,jr)           
      yr=verty(ir,jr)
      xr1=vertx(ir1,jr1)
      yr1=verty(ir1,jr1)

      ur=Vstarx(ir,jr)    ! u star r
      vr=Vstary(ir,jr)    ! v star r
      ur1=Vstarx(ir1,jr1)
      vr1=Vstary(ir1,jr1)

      dAdt=dAdt+(ur*yr1+vr1*xr-ur1*yr-vr*xr1)*0.5
      !c3
      ir=i
      jr=j
      ir1=i-1
      jr1=j

      xr=vertx(ir,jr)           
      yr=verty(ir,jr)
      xr1=vertx(ir1,jr1)
      yr1=verty(ir1,jr1)

      ur=Vstarx(ir,jr)    ! u star r
      vr=Vstary(ir,jr)    ! v star r
      ur1=Vstarx(ir1,jr1)
      vr1=Vstary(ir1,jr1)

      dAdt=dAdt+(ur*yr1+vr1*xr-ur1*yr-vr*xr1)*0.5
      !c4
      ir=i-1
      jr=j
      ir1=i-1
      jr1=j-1

      xr=vertx(ir,jr)           
      yr=verty(ir,jr)
      xr1=vertx(ir1,jr1)
      yr1=verty(ir1,jr1)

      ur=Vstarx(ir,jr)    ! u star r
      vr=Vstary(ir,jr)    ! v star r
      ur1=Vstarx(ir1,jr1)
      vr1=Vstary(ir1,jr1)

      dAdt=dAdt+(ur*yr1+vr1*xr-ur1*yr-vr*xr1)*0.5

      tV=min(CV*amass(i,j)*tau(i,j)/abs(dAdt),tV)   

   enddo
   enddo
   
   dt= min(tV,tE) * cfl
    
 end
