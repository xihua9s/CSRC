
!> 20231130
!> Xihua, xu_xihua@iapcm.ac.cn
!> A parameter-free staggered-grid Lagragian scheme for two-dimensional
!> compressible flow problems 
!> Criticism is welcome.

Program main

   call mjsls_ca()
   
end program


subroutine mjsls_ca()
	
	implicit none
	real*8 ::gamar
	integer :: nx,ny
	parameter (nx=70,ny=30)

	real*8  :: den(1:nx,1:ny),pre(1:nx,1:ny),amass(1:nx,1:ny)
	real*8  :: sound(1:nx,1:ny),enin(1:nx,1:ny)	
	real*8  :: vertx(0:nx,0:ny),verty(0:nx,0:ny)
	real*8  :: Vstarx(0:nx,0:ny),Vstary(0:nx,0:ny)
	real*8  :: Vertex_Type_norm_pre(0:nx,0:ny,3)
	integer :: eos(1:nx,1:ny),bnd_type(0:nx,0:ny)
   
   
	real*8  :: supt,dt,t,cfl,g
	integer :: it,t_out,it_out

	g = 0.0d0   !> gravity

	call initial_sls_triple(supt,den,pre,sound,enin,eos,&
			vertx,verty,amass,gamar,nx,ny,vstarx,vstary,bnd_type,Vertex_Type_norm_pre)         

	t      = 0.d0 
	it     = 0
   cfl    = 0.2d0
   t_out  = 10
   it_out = 1
   
	do while(t.lt.supt)

		call getdt0(sound,amass,den,vertx,verty,Vstarx,Vstary,nx,ny,cfl,dt)
      
		call cal_mjsls_details(nx,ny,den,pre,sound,enin,eos,vertx,verty,amass, &
           gamar,vstarx,vstary,dt,bnd_type,Vertex_Type_norm_pre,g)

		if(t+dt.gt.supt) dt=supt-t
		t  = t  + dt
		it = it + 1
		write(*,*) "it=",it,"T=",t,"dt=",dt 
      
      if( t > supt/t_out*it_out ) then  
			call output_sls(nx,ny,den,pre,enin,vertx,verty,gamar,vstarx,vstary,it_out) 
         it_out = it_out + 1
		endif 

	enddo       

	call output_sls_final(nx,ny,den,pre,enin,vertx,verty,gamar,vstarx,vstary,it) 

end subroutine mjsls_ca
