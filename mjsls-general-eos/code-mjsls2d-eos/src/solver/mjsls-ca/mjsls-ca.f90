
!> Minor Jump Staggered Lagrangian Scheme in Cylindral

module mdu_mjsls_ca

   use mdu_update
contains 

subroutine mjsls_ca(ini_cmt,nx,ny)
	use common_inicmt
	implicit none
	
   type(inicmt) :: ini_cmt
	integer :: nx,ny
	
	real*8  :: den(1:nx,1:ny),pre(1:nx,1:ny),amass(1:nx,1:ny)
	real*8  :: sound(1:nx,1:ny),enin(1:nx,1:ny),smass(1:nx,1:ny,1:4)
	real*8  :: smass_tri(1:nx,1:ny,1:8)
	real*8  :: vertx(0:nx,0:ny),verty(0:nx,0:ny),pmass(0:nx,0:ny)
	real*8  :: Vstarx(0:nx,0:ny),Vstary(0:nx,0:ny)
	real*8  :: Vertex_Type_norm_pre(0:nx,0:ny,3)
	integer :: eos(1:nx,1:ny),bnd_type(0:nx,0:ny)
	real*8 ::gamar
      
	real*8  :: supt,dt,t,cfl,g,area
	integer :: it,t_out,it_out
	real*8  :: l1norm,dx,dy,ddx,ddy,l2norm,l8norm
	integer :: i,j,n
	character*99	:: case_name,nodal_solver

	t      = ini_cmt%it 
   cfl    = ini_cmt%cfl
   t_out  = ini_cmt%t_out
	case_name = ini_cmt%case_name 
   nodal_solver = ini_cmt%node_solver
	
   it     = 0
	supt   = 0.d0 
   it_out = 1

	call initial_ca(supt,den,pre,sound,enin,eos, &
   vertx,verty,amass,gamar,nx,ny,vstarx,vstary,bnd_type,&
   Vertex_Type_norm_pre,pmass,smass,g,case_name,smass_tri)

	call output_init(nx,ny,den,pre,enin,vertx,verty,gamar,vstarx,vstary,it) 
   
	do while(t.lt.supt)

		call getdt0(sound,amass,den,vertx,verty,Vstarx,Vstary,nx,ny,cfl,dt,it)
      
      select case (nodal_solver)

      case ('prine')

   		call cal_mjsls_details_ca(nx,ny,den,pre,sound,enin,eos,vertx,verty,amass, &
              gamar,vstarx,vstary,dt,bnd_type,Vertex_Type_norm_pre,g,pmass,smass,smass_tri)

      case ('prine_dt')

   		call cal_mjsls_details_ca_dpdtaudt(nx,ny,den,pre,sound,enin,eos,vertx,verty,amass, &
              gamar,vstarx,vstary,dt,bnd_type,Vertex_Type_norm_pre,g,pmass,smass,smass_tri)

      case default

         print*, 'nodal_solver=', nodal_solver 

      end select

		if(t+dt.gt.supt) dt=supt-t
		t  = t  + dt
		it = it + 1
		write(*,*) "it=",it,"T=",t,"dt=",dt 

      !if( t > supt/t_out*it_out .or. mod(ini_cmt%t_out,it) == 0) then  
      if( t > supt/t_out*it_out) then  
			call output_it_time(nx,ny,den,pre,enin,vertx,verty,gamar,vstarx,vstary,it_out) 
         it_out = it_out + 1
		endif 

      if( mod(it,ini_cmt%it_out) == 0) then  
			call output_it_time(nx,ny,den,pre,enin,vertx,verty,gamar,vstarx,vstary,it) 
		endif       

	enddo       

	call output_final(nx,ny,den,pre,enin,vertx,verty,gamar,vstarx,vstary,it) 

end subroutine mjsls_ca



!> 20230114
!> Xihua 
!> cal_mjsls_details_cy 
subroutine cal_mjsls_details_ca(nx,ny,den,pre,sound,enin,eos,vertx,verty,amass, &
            gamar,vstarx,vstary,dt,bnd_type,Vertex_Type_norm_pre,g,pmass,smass,smass_tri)
   implicit none 
   
	integer,intent(in) :: nx,ny
	real*8 :: den(1:nx,1:ny),pre(1:nx,1:ny),amass(1:nx,1:ny)
	real*8 :: sound(1:nx,1:ny),enin(1:nx,1:ny)
	real*8 :: xc(1:nx,1:ny),yc(1:nx,1:ny),g,pmass(0:nx,0:ny)
   real*8 :: pre1(4,1:nx,1:ny),pre1_in(4,1:nx,1:ny)
	integer :: eos(1:nx,1:ny),bnd_type(0:nx,0:ny)
	real*8  :: Vertex_Type_norm_pre(0:nx,0:ny,3),smass(1:nx,1:ny,1:4)
   real*8  :: smass_tri(1:nx,1:ny,1:8)
   
	real*8 :: vertx(0:nx,0:ny),verty(0:nx,0:ny)
	real*8 :: vstarx(0:nx,0:ny),vstary(0:nx,0:ny)
	real*8 :: velx(0:nx,0:ny),vely(0:nx,0:ny)
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
   call cal_mjsls_details_subcell_pre4_in_ca(nx,ny,eos,vertx, &
         verty,vstarx,vstary,den,pre,sound,gamar,amass,xc,yc, &
         pre1,pre1_in,smass,smass_tri)
                  
   !> calculate the force of each nodal
   call cal_node_force_ca(nx,ny,amass,pre1,pre1_in,vertx,verty, &
         xc,yc,bnd_type,Vertex_Type_norm_pre,pmass,fx_pre,fy_pre)      

   !> calculate the velocity of each nodal
   call cal_node_velocity_ca(nx,ny,pmass,dt,bnd_type, &
         Vertex_Type_norm_pre,vstarx,vstary,vstarx_new,&
         vstary_new,fx_pre,fy_pre,g,vertx,verty,smass)         

   !> calculate the average velocity and update the vertex 
   velx = (vstarx + vstarx_new) / 2.d0 
   vely = (vstary + vstary_new) / 2.d0  

   call update_vertex(nx,ny,vertx,verty,vertx_new,verty_new,velx,vely,dt)

   !> update density and internal energy
   call update_den_enin_ca(nx,ny,eos,vertx,verty,&
         vertx_new,verty_new,velx,vely,dt,pre1,den, &
         enin,amass,gamar,den_new,enin_new,pre_new,sound_new,g) 

   ! vertx  = vertx_new
   ! verty  = verty_new 
   ! vstarx = vstarx_new
   ! vstary = vstary_new       
   ! den    = den_new
   ! pre    = pre_new
   ! sound  = sound_new
   ! enin   = enin_new 
   ! return 
   !> cor step -------------------------------
    
   !> calculate pre2 and pre2_in of each subcell 
   call cal_mjsls_details_subcell_pre4_in_ca(nx,ny,eos,&
         vertx_new,verty_new,vstarx_new,vstary_new, &
         den_new,pre_new,sound_new,gamar,amass, &
         xc_pre,yc_pre,pre2,pre2_in,smass,smass_tri)

   !> calculate the force of each nodal 
   call cal_node_force_ca(nx,ny,amass,(pre1+pre2)/2.d0, &
         (pre1_in+pre2_in)/2.d0,vertx_new,verty_new, &
         xc_pre,yc_pre,bnd_type,Vertex_Type_norm_pre,&
         pmass,fx_cor,fy_cor)         

   !> calculate the velocity of each nodal
   call cal_node_velocity_ca(nx,ny,pmass,dt,bnd_type, &
         Vertex_Type_norm_pre,vstarx,vstary,vstarx_new, &
         vstary_new,(fx_pre+fx_cor)/2.d0,(fy_pre+fy_cor)/2.d0,g,vertx,verty,smass)                

   !> calculate the average velocity and update the vertex 
   velx = (vstarx + vstarx_new) / 2.d0 
   vely = (vstary + vstary_new) / 2.d0   
   call update_vertex(nx,ny,vertx,verty,vertx_new,verty_new,velx,vely,dt)

   !> update density and internal energy
   call update_den_enin_ca(nx,ny,eos,vertx,verty,vertx_new,verty_new,velx,vely,dt,&
            (pre1+pre2)/2.d0,den,enin,amass,gamar,den_new,enin_new,pre_new,sound_new,g)    
        
   vertx  = vertx_new
   verty  = verty_new 
   vstarx = vstarx_new
   vstary = vstary_new       
   den    = den_new
   pre    = pre_new
   sound  = sound_new
   enin   = enin_new    

end subroutine cal_mjsls_details_ca


!> 20230114
!> Xihua 
!> cal_mjsls_details_cy 
subroutine cal_mjsls_details_ca_dpdtaudt(nx,ny,den,pre,sound,enin,eos,vertx,verty,amass, &
            gamar,vstarx,vstary,dt,bnd_type,Vertex_Type_norm_pre,g,pmass,smass,smass_tri)
   implicit none 
   
	integer,intent(in) :: nx,ny
	real*8 :: den(1:nx,1:ny),pre(1:nx,1:ny),amass(1:nx,1:ny)
	real*8 :: sound(1:nx,1:ny),enin(1:nx,1:ny)
	real*8 :: xc(1:nx,1:ny),yc(1:nx,1:ny),g,pmass(0:nx,0:ny)
   real*8 :: pre1(4,1:nx,1:ny),pre1_in(4,1:nx,1:ny)
	integer :: eos(1:nx,1:ny),bnd_type(0:nx,0:ny)
	real*8  :: Vertex_Type_norm_pre(0:nx,0:ny,3),smass(1:nx,1:ny,1:4)
   real*8  :: smass_tri(1:nx,1:ny,1:8)
   
	real*8 :: vertx(0:nx,0:ny),verty(0:nx,0:ny)
	real*8 :: vstarx(0:nx,0:ny),vstary(0:nx,0:ny)
	real*8 :: velx(0:nx,0:ny),vely(0:nx,0:ny)
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
   call cal_mjsls_details_subcell_pre4_in_ca_dpdtaudt(nx,ny,eos,vertx, &
         verty,vstarx,vstary,den,pre,sound,gamar,amass,xc,yc, &
         pre1,pre1_in,smass,smass_tri,dt)
                  
   !> calculate the force of each nodal
   call cal_node_force_ca(nx,ny,amass,pre1,pre1_in,vertx,verty, &
         xc,yc,bnd_type,Vertex_Type_norm_pre,pmass,fx_pre,fy_pre)      

   !> calculate the velocity of each nodal
   call cal_node_velocity_ca(nx,ny,pmass,dt,bnd_type, &
         Vertex_Type_norm_pre,vstarx,vstary,vstarx_new,&
         vstary_new,fx_pre,fy_pre,g,vertx,verty,smass)         

   !> calculate the average velocity and update the vertex 
   velx = (vstarx + vstarx_new) / 2.d0 
   vely = (vstary + vstary_new) / 2.d0  

   call update_vertex(nx,ny,vertx,verty,vertx_new,verty_new,velx,vely,dt)

   !> update density and internal energy
   call update_den_enin_ca(nx,ny,eos,vertx,verty,&
         vertx_new,verty_new,velx,vely,dt,pre1,den, &
         enin,amass,gamar,den_new,enin_new,pre_new,sound_new,g) 

   ! vertx  = vertx_new
   ! verty  = verty_new 
   ! vstarx = vstarx_new
   ! vstary = vstary_new       
   ! den    = den_new
   ! pre    = pre_new
   ! sound  = sound_new
   ! enin   = enin_new 
   ! return 
   !> cor step -------------------------------
    
   !> calculate pre2 and pre2_in of each subcell 
   call cal_mjsls_details_subcell_pre4_in_ca_dpdtaudt(nx,ny,eos,&
         vertx_new,verty_new,vstarx_new,vstary_new, &
         den_new,pre_new,sound_new,gamar,amass, &
         xc_pre,yc_pre,pre2,pre2_in,smass,smass_tri,dt)

   !> calculate the force of each nodal 
   call cal_node_force_ca(nx,ny,amass,(pre1+pre2)/2.d0, &
         (pre1_in+pre2_in)/2.d0,vertx_new,verty_new, &
         xc_pre,yc_pre,bnd_type,Vertex_Type_norm_pre,&
         pmass,fx_cor,fy_cor)         

   !> calculate the velocity of each nodal
   call cal_node_velocity_ca(nx,ny,pmass,dt,bnd_type, &
         Vertex_Type_norm_pre,vstarx,vstary,vstarx_new, &
         vstary_new,(fx_pre+fx_cor)/2.d0,(fy_pre+fy_cor)/2.d0,g,vertx,verty,smass)                

   !> calculate the average velocity and update the vertex 
   velx = (vstarx + vstarx_new) / 2.d0 
   vely = (vstary + vstary_new) / 2.d0   
   call update_vertex(nx,ny,vertx,verty,vertx_new,verty_new,velx,vely,dt)

   !> update density and internal energy
   call update_den_enin_ca(nx,ny,eos,vertx,verty,vertx_new,verty_new,velx,vely,dt,&
            (pre1+pre2)/2.d0,den,enin,amass,gamar,den_new,enin_new,pre_new,sound_new,g)    
        
   vertx  = vertx_new
   verty  = verty_new 
   vstarx = vstarx_new
   vstary = vstary_new       
   den    = den_new
   pre    = pre_new
   sound  = sound_new
   enin   = enin_new    

end subroutine cal_mjsls_details_ca_dpdtaudt



end module mdu_mjsls_ca