
module mdu_mjsls_cy 
   use mdu_update

   implicit none 
contains 

!> Minor Jump Staggered Lagrangian Scheme in Cylindral
subroutine mjsls_cy_aw()
	
	implicit none
	real*8 ::gamar
	integer :: nx,ny
	parameter (nx=100,ny=10)

	real*8  :: den(1:nx,1:ny),pre(1:nx,1:ny),amass(1:nx,1:ny)
	real*8  :: sound(1:nx,1:ny),enin(1:nx,1:ny),smass(1:nx,1:ny,1:4),smass_tri(1:nx,1:ny,1:8)
	real*8  :: vertx(0:nx,0:ny),verty(0:nx,0:ny),pmass(0:nx,0:ny)
	real*8  :: Vstarx(0:nx,0:ny),Vstary(0:nx,0:ny)
	real*8  :: Vertex_Type_norm_pre(0:nx,0:ny,3)
	integer :: eos(1:nx,1:ny),bnd_type(0:nx,0:ny)
      
	real*8  :: supt,dt,t,cfl,g,area
	integer :: it,t_out,it_out
	real*8  :: l1norm,dx,dy,ddx,ddy,l2norm,l8norm
	integer :: i,j,n

   supt = 0.d0 
	! call initial_ca(supt,den,pre,sound,enin,eos, &
   ! vertx,verty,amass,gamar,nx,ny,vstarx,vstary,bnd_type,&
   ! Vertex_Type_norm_pre,pmass,smass,g)

	call output_init(nx,ny,den,pre,enin,vertx,verty,gamar,vstarx,vstary,it) 

	t      = 0.d0 
	it     = 0
   cfl    = 0.2d0
   t_out  = 10
   it_out = 1
   
	do while(t.lt.supt)

		call getdt0(sound,amass,den,vertx,verty,Vstarx,Vstary,nx,ny,cfl,dt,it)
      
		call mjsls_cy_aw_details(nx,ny,den,pre,sound,enin,eos,vertx,verty,amass, &
           gamar,vstarx,vstary,dt,bnd_type,Vertex_Type_norm_pre,g,pmass,smass)

		if(t+dt.gt.supt) dt=supt-t
		t  = t  + dt
		it = it + 1
		write(*,*) "it=",it,"T=",t,"dt=",dt 
      
      if( t > supt/t_out*it_out ) then  
			call output_it_time(nx,ny,den,pre,enin,vertx,verty,gamar,vstarx,vstary,it_out) 
         it_out = it_out + 1
		endif 

	enddo       

	call output_final(nx,ny,den,pre,enin,vertx,verty,gamar,vstarx,vstary,it) 

end subroutine mjsls_cy_aw


!> 20230114
!> Xihua 
!> mjsls_cy_aw_details 
subroutine mjsls_cy_aw_details(nx,ny,den,pre,sound,enin,eos,vertx,verty,amass, &
            gamar,vstarx,vstary,dt,bnd_type,Vertex_Type_norm_pre,g,pmass,smass)
   implicit none 
   
	integer,intent(in) :: nx,ny
	real*8 :: den(1:nx,1:ny),pre(1:nx,1:ny),amass(1:nx,1:ny)
	real*8 :: sound(1:nx,1:ny),enin(1:nx,1:ny)
	real*8 :: xc(1:nx,1:ny),yc(1:nx,1:ny),g,pmass(0:nx,0:ny)
   real*8 :: pre1(4,1:nx,1:ny),pre1_in(4,1:nx,1:ny)
	integer :: eos(1:nx,1:ny),bnd_type(0:nx,0:ny)
	real*8  :: Vertex_Type_norm_pre(0:nx,0:ny,3),smass(1:nx,1:ny,1:4),smass_tri(1:nx,1:ny,1:8)
   
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

   !> pre step -------------------------------
   
   !> calculate pre1 and pre1_in of each subcell 
   call cal_mjsls_details_subcell_pre4_in_cy(nx,ny,eos,vertx, &
         verty,vstarx,vstary,den,pre,sound,gamar,amass,xc,yc, &
         pre1,pre1_in)
            
   !> calculate the force of each nodal
   call cal_node_force_cy(nx,ny,amass,pre1,pre1_in,vertx,verty, &
         xc,yc,bnd_type,Vertex_Type_norm_pre,mass_p,fx_pre,fy_pre)      

   !> calculate the velocity of each nodal
   call cal_node_velocity_cy(nx,ny,pmass,dt,bnd_type, &
         Vertex_Type_norm_pre,vstarx,vstary,vstarx_new,&
         vstary_new,fx_pre,fy_pre,g,vertx,verty,smass)         

   ! print*, vstarx(:,0) ; pause 29
   ! print*, vstarx_new(:,0) ; pause 60
   ! print*, vstarx(:,1) ; pause 28
   ! print*, vstarx_new(:,1) ; pause 23

   !> calculate the average velocity and update the vertex 
   velx = (vstarx + vstarx_new) / 2.d0 
   vely = (vstary + vstary_new) / 2.d0   
   call update_vertex(nx,ny,vertx,verty,vertx_new,verty_new,velx,vely,dt)

   !> update density and internal energy
   call update_den_enin_cy(nx,ny,eos,vertx,verty,&
         vertx_new,verty_new,velx,vely,dt,pre1,den, &
         enin,amass,gamar,den_new,enin_new,pre_new,sound_new,g) 

   vertx  = vertx_new
   verty  = verty_new 
   vstarx = vstarx_new
   vstary = vstary_new       
   den    = den_new
   pre    = pre_new
   sound  = sound_new
   enin   = enin_new 
   return 
   !> cor step -------------------------------
    
   !> calculate pre2 and pre2_in of each subcell 
   call cal_mjsls_details_subcell_pre4_in_cy(nx,ny,eos,&
         vertx_new,verty_new,vstarx_new,vstary_new, &
         den_new,pre_new,sound_new,gamar,amass, &
         xc_pre,yc_pre,pre2,pre2_in)

   !> calculate the force of each nodal 
   call cal_node_force_cy(nx,ny,amass,(pre1+pre2)/2.d0, &
         (pre1_in+pre2_in)/2.d0,vertx_new,verty_new, &
         xc_pre,yc_pre,bnd_type,Vertex_Type_norm_pre,&
         mass_p,fx_cor,fy_cor)         

   !> calculate the velocity of each nodal
   call cal_node_velocity_cy(nx,ny,mass_p,dt,bnd_type, &
         Vertex_Type_norm_pre,vstarx,vstary,vstarx_new, &
         vstary_new,(fx_pre+fx_cor)/2.d0,(fy_pre+fy_cor)/2.d0,g,vertx,verty,smass)                

   !> calculate the average velocity and update the vertex 
   velx = (vstarx + vstarx_new) / 2.d0 
   vely = (vstary + vstary_new) / 2.d0   
   call update_vertex(nx,ny,vertx,verty,vertx_new,verty_new,velx,vely,dt)

   !> update density and internal energy
   call update_den_enin_cy(nx,ny,eos,vertx,verty,vertx_new,verty_new,velx,vely,dt,&
            (pre1+pre2)/2.d0,den,enin,amass,gamar,den_new,enin_new,pre_new,sound_new,g)    
        
   vertx  = vertx_new
   verty  = verty_new 
   vstarx = vstarx_new
   vstary = vstary_new       
   den    = den_new
   pre    = pre_new
   sound  = sound_new
   enin   = enin_new    

end subroutine mjsls_cy_aw_details

end module mdu_mjsls_cy