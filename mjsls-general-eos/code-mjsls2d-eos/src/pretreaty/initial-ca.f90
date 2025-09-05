
subroutine initial_ca(supt,den,pre,sound,enin,eos, &
   vertx,verty,amass,gamar,nx,ny,vstarx,vstary,bnd_type,&
   Vertex_Type_norm_pre,pmass,smass,g,case_name,smass_tri)
	implicit none
	real*8 :: den(1:nx,1:ny),pre(1:nx,1:ny),amass(1:nx,1:ny)
	real*8 :: sound(1:nx,1:ny),enin(1:nx,1:ny)
	real*8 :: area(1:nx,1:ny),volm(1:nx,1:ny)
	real*8 :: vertx(0:nx,0:ny),verty(0:nx,0:ny),pmass(0:nx,0:ny)
	real*8 :: vstarx(0:nx,0:ny),vstary(0:nx,0:ny)
	integer  :: eos(1:nx,1:ny),bnd_type(0:nx,0:ny)
	real*8  :: Vertex_Type_norm_pre(0:nx,0:ny,3)
	integer :: nx,ny
	real*8  :: gamar,supt,ca_area,cy_volm,parea(0:nx,0:ny),pvolm(0:nx,0:ny)
   real*8  :: subcell_volm_cy(1:nx,1:ny,1:4),smass(1:nx,1:ny,1:4)
   real*8  :: smass_tri(1:nx,1:ny,1:8)
   real*8  :: subcell_area_ca(1:nx,1:ny,1:4) ,g
   real*8  :: subcell_area_cy_aw(1:nx,1:ny,1:4)    
   real*8  :: subcell_volm_cy_aw(1:nx,1:ny,1:4) 
	character*99	:: case_name   

	g = 0.0d0   !> gravity

   select case(case_name)
   case('sod_rec')
      call initial_sod_rectangle_ca(supt,den,pre,sound,enin,eos, &
            vertx,verty,amass,gamar,nx,ny,vstarx,vstary, &
            bnd_type,Vertex_Type_norm_pre,pmass,smass,smass_tri)
   case('s123_rec')
      call initial_s123_rectangle_ca(supt,den,pre,sound,enin,eos, &
            vertx,verty,amass,gamar,nx,ny,vstarx,vstary, &
            bnd_type,Vertex_Type_norm_pre,pmass,smass,smass_tri)
   case('exp2vac_rec')
      call initial_exp2vac_rectangle_ca(supt,den,pre,sound,enin,eos, &
            vertx,verty,amass,gamar,nx,ny,vstarx,vstary, &
            bnd_type,Vertex_Type_norm_pre,pmass,smass,smass_tri)
   case('exp2vac_cir')
      call initial_exp2vac_circle_ca(supt,den,pre,sound,enin,eos, &
            vertx,verty,amass,gamar,nx,ny,vstarx,vstary, &
            bnd_type,Vertex_Type_norm_pre,pmass,smass,smass_tri)
   case('leblanc_rec')
      call initial_leblanc_rectangle_ca(supt,den,pre,sound,enin,eos, &
            vertx,verty,amass,gamar,nx,ny,vstarx,vstary, &
            bnd_type,Vertex_Type_norm_pre,pmass,smass,smass_tri)
   case('sedov_rec')
      call initial_sedov_rectangle_ca(supt,den,pre,sound,enin,eos, &
            vertx,verty,amass,gamar,nx,ny,vstarx,vstary, &
            bnd_type,Vertex_Type_norm_pre,pmass,smass,smass_tri)
   case('sedov4_rec')
      call initial_sedov4_rectangle_ca(supt,den,pre,sound,enin,eos, &
            vertx,verty,amass,gamar,nx,ny,vstarx,vstary, &
            bnd_type,Vertex_Type_norm_pre,pmass,smass,smass_tri)
   case('noh_rec')
      call initial_noh_rectangle_ca(supt,den,pre,sound,enin,eos, &
            vertx,verty,amass,gamar,nx,ny,vstarx,vstary, &
            bnd_type,Vertex_Type_norm_pre,pmass,smass,smass_tri)
   case('noh_cir')
      call initial_noh_circle_ca(supt,den,pre,sound,enin,eos, &
            vertx,verty,amass,gamar,nx,ny,vstarx,vstary, &
            bnd_type,Vertex_Type_norm_pre,pmass,smass,smass_tri)
   case('saltzmann')
      call initial_saltzmann_ca(supt,den,pre,sound,enin,eos, &
            vertx,verty,amass,gamar,nx,ny,vstarx,vstary, &
            bnd_type,Vertex_Type_norm_pre,pmass,smass,smass_tri)       
   case('triple')
      call initial_triple_ca(supt,den,pre,sound,enin,eos, &
            vertx,verty,amass,gamar,nx,ny,vstarx,vstary, &
            bnd_type,Vertex_Type_norm_pre,pmass,smass,smass_tri)    
   case('rti')
      call initial_rti_ca(supt,den,pre,sound,enin,eos, &
            vertx,verty,amass,gamar,nx,ny,vstarx,vstary, &
            bnd_type,Vertex_Type_norm_pre,pmass,smass,g)  
   case('acousticpulse')
      call initial_acousticpulse_ca(supt,den,pre,sound,enin,eos, &
            vertx,verty,amass,gamar,nx,ny,vstarx,vstary, &
            bnd_type,Vertex_Type_norm_pre,pmass,smass,smass_tri)                
   case  default
      print*, 'case_name:',case_name
      stop 'case_name is error'
   end select
end subroutine initial_ca



