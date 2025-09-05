!> Last Modified : Sat Feb 25 16:08:50 2023
module mdu_fluid_euler_eqn_solver_centered_update

	use common_inicmt
	use common_meshphys
	use fluid_euler_eqn_basic
	use MDU_Initial
	use dt
	use output
	implicit none

contains


	Subroutine update_vertex(ini_cmt,mesh_phys)
		implicit none

		type(inicmt)   :: ini_cmt 
		type(meshphys) :: mesh_phys 

		integer :: ix,nu

		select case(ini_cmt%cacysp)
		case('Cartesian')
			nu = 0
		case('Cylindrical')
			nu = 1
		case('Spheral')
			nu = 2
		case default
			print*, 'ini_cmt%cacysp=', ini_cmt%cacysp
			stop 'ini_cmt%cacysp is error.'
		end select      

		do ix = 1 , ini_cmt%nx
		mesh_phys%vol_size_old(ix) = mesh_phys%vertex(ix+1)**(nu+1)/(nu+1) - &
			mesh_phys%vertex(ix)**(nu+1)/(nu+1)
		enddo    

		Do ix = 1 , ini_cmt%nx+1
		mesh_phys%vertex(ix) = mesh_phys%vertex(ix) + mesh_phys%vel_p(ix) * ini_cmt%dt
		Enddo
		Do ix = 1 , ini_cmt%nx
		mesh_phys%vol_size_new(ix) = mesh_phys%vertex(ix+1)**(nu+1)/(nu+1) - &
			mesh_phys%vertex(ix)**(nu+1)/(nu+1)         
		Enddo   

	End Subroutine update_vertex   

	Subroutine update_phys(ini_cmt,mesh_phys)
		implicit none

		type(inicmt)   :: ini_cmt 
		type(meshphys) :: mesh_phys

		Select Case (ini_cmt%node_solver)
		Case('Acoustic')

			call update_phys_details(ini_cmt,mesh_phys)

		Case('MJCLS')

			call update_phys_details(ini_cmt,mesh_phys)      

		Case('MJCLS2P2H')

			call update_phys_details(ini_cmt,mesh_phys)      

		Case('UeAS')

			call update_phys_Xihua_dudedE(ini_cmt,mesh_phys)

		Case default

			print*, 'ini_cmt%node_solver is', ini_cmt%node_solver
			stop 'ini_cmt%node_solver is Error, pls Call Xihua.'

		End Select

	End Subroutine update_phys


	Subroutine update_phys_details(ini_cmt,mesh_phys)
		implicit none

		type(inicmt)   :: ini_cmt 
		type(meshphys) :: mesh_phys

		integer :: ix,nu
		real*8  :: den,vel,ene,cvar(3),pvar(3),sos,ein,gamma
		real*8  :: r1,r2,rbar

		select case(ini_cmt%cacysp)
		case('Cartesian')
			nu = 0
		case('Cylindrical')
			nu = 1
		case('Spheral')
			nu = 2
		case default
			print*, 'ini_cmt%cacysp=', ini_cmt%cacysp
			stop 'ini_cmt%cacysp is error.'
		end select

		Do ix = 1 , ini_cmt%nx
		den = mesh_phys%con_phys(ix,1)
		vel = mesh_phys%con_phys(ix,2) / mesh_phys%con_phys(ix,1)
		ene = mesh_phys%con_phys(ix,3) / mesh_phys%con_phys(ix,1)

		r1 = mesh_phys%vertex(ix)
		r2 = mesh_phys%vertex(ix+1)
		rbar = ( r2**(nu+2)/(nu+2) - r1**(nu+2)/(nu+2) ) / ( r2**(nu+1)/(nu+1)-r1**(nu+1)/(nu+1) )

		mesh_phys%con_phys(ix,1) = mesh_phys%mass(ix) / mesh_phys%vol_size_new(ix)
		vel = vel - rbar**(nu)*(mesh_phys%pstar(ix+1)-mesh_phys%pstar(ix))*ini_cmt%dt/mesh_phys%mass(ix)
		ene = ene - (mesh_phys%pstar(ix+1)*mesh_phys%vel_p(ix+1)*r2**(nu) - &
			mesh_phys%pstar(ix)  *mesh_phys%vel_p(ix)  *r1**(nu))*ini_cmt%dt / mesh_phys%mass(ix)
		mesh_phys%con_phys(ix,2) = mesh_phys%con_phys(ix,1) * vel
		mesh_phys%con_phys(ix,3) = mesh_phys%con_phys(ix,1) * Ene

		cvar = mesh_phys%con_phys(ix,1:3)

		Call ieos_to_gamma(mesh_phys%eos(ix),gamma)
		Call CVarToPVar(3, CVar,PVar,gamma,mesh_phys%eos(ix))     
		Call PVarToOther(3, PVar,sos,ein,gamma,mesh_phys%eos(ix))

		mesh_phys%vol_size_old(ix) = mesh_phys%vertex(ix+1)**(nu+1)/(nu+1) - &
			mesh_phys%vertex(ix)**(nu+1)/(nu+1)
		mesh_phys%pre(ix) = pvar(3)
		mesh_phys%sos(ix) = sos
		mesh_phys%ein(ix) = ein

		enddo


	End Subroutine update_phys_details


	Subroutine update_phys_Xihua_dudedE(ini_cmt,mesh_phys)
		implicit none

		type(inicmt)   :: ini_cmt 
		type(meshphys) :: mesh_phys

		integer :: ix,nn
		real*8  :: den,vel,ene,cvar(3),pvar(3),sos,ein,gamma,vv
		real*8  :: q , b , tmp_pre , am , veltmp1, veltmp2, veltmp
		real*8  :: pstar_L,ustar_L,dstar_L,sstar_L
		real*8  :: pstar_R,ustar_R,dstar_R,sstar_R,pre_bar_L,pre_bar_R
		real*8  :: delta_e,delta_ptau,epsi,v0,v1,v2,alpha,beta,r1,r2,rbar
		integer :: nu,signa
		real*8  :: sigma

		select case(ini_cmt%cacysp)
		case('Cartesian')
			nu = 0
		case('Cylindrical')
			nu = 1
		case('Spheral')
			nu = 2
		case default
			print*, 'ini_cmt%cacysp=', ini_cmt%cacysp
			stop 'ini_cmt%cacysp is error.'
		end select   

		!> nn = 1000000000
		sigma = 1.d-8
		Do ix = 1 , ini_cmt%nx
		signa = 0
		!> n ʱ�̵�ֵ
		den = mesh_phys%con_phys(ix,1)
		vel = mesh_phys%con_phys(ix,2) / mesh_phys%con_phys(ix,1)
		ene = mesh_phys%con_phys(ix,3) / mesh_phys%con_phys(ix,1)
		ein = ene - vel*vel/2.d0

		r1 = mesh_phys%vertex(ix)
		r2 = mesh_phys%vertex(ix+1)
		rbar = ( r2**(nu+2)/(nu+2) - r1**(nu+2)/(nu+2) ) / ( r2**(nu+1)/(nu+1)-r1**(nu+1)/(nu+1) )

		!> n + 1 ʱ�̵�ֵ
		mesh_phys%con_phys(ix,1) = mesh_phys%mass(ix) / mesh_phys%vol_size_new(ix)


		! veltmp  = (nn/2.d0-1.d0)/(nn/2.d0) * vel + &
		! 1.d0/nn * mesh_phys%vel_p(ix+1) + 1.d0/nn * mesh_phys%vel_p(ix+0) 
		veltmp  = (1.d0 - 2.d0*sigma) * vel + &
			sigma * mesh_phys%vel_p(ix+1) + sigma * mesh_phys%vel_p(ix+0)


		vel = vel - rbar**(nu)*(mesh_phys%pstar(ix+1)-mesh_phys%pstar(ix))*ini_cmt%dt/mesh_phys%mass(ix)
		ene = ene - (mesh_phys%pstar(ix+1)*mesh_phys%vel_p(ix+1)*r2**(nu) - &
			mesh_phys%pstar(ix)  *mesh_phys%vel_p(ix)  *r1**(nu))*ini_cmt%dt / mesh_phys%mass(ix)      
		ein = ein - (  mesh_phys%pstar(ix+1) * mesh_phys%vel_p(ix+1) *r2**(nu) - &
			mesh_phys%pstar(ix+1) * veltmp * rbar**(nu)             - &
			mesh_phys%pstar(ix+0) * mesh_phys%vel_p(ix+0) *r1**(nu) + &
			mesh_phys%pstar(ix+0) * veltmp * rbar**(nu)                 )  * &
			ini_cmt%dt / mesh_phys%mass(ix)      


		mesh_phys%ein(ix) = ein

		vv = (ene - ein) * 2.d0   
		!> ȡ����
		if(vv < 0.d0) then  
			print*, vv , ix
			signa = 1
			!write(187,*),ix,vv,ini_cmt%current_time
			vv = -vv
		else
			vv = sqrt(vv)
		endif

		!> ȡ vel �ķ���
		if( vel > 0.d0) then
			vel = vv
		else
			vel = -vv
		endif         

		mesh_phys%con_phys(ix,2) = mesh_phys%con_phys(ix,1) * vel
		mesh_phys%con_phys(ix,3) = mesh_phys%con_phys(ix,1) * Ene

		if(signa == 1) then  !> �����E-e<0ʱ����Ԫ�ٶȷ����Ƿ���
			!write(188,*), vel,mesh_phys%con_phys(ix,2)
		endif

		cvar = mesh_phys%con_phys(ix,1:3)

		Call ieos_to_gamma(mesh_phys%eos(ix),gamma)
		Call CVarToPVar(3, CVar,PVar,gamma,mesh_phys%eos(ix))     
		Call PVarToOther(3, PVar,sos,ein,gamma,mesh_phys%eos(ix))

		mesh_phys%vol_size_old(ix) = mesh_phys%vertex(ix+1)**(nu+1)/(nu+1) - &
			mesh_phys%vertex(ix)**(nu+1)/(nu+1)
		mesh_phys%pre(ix) = pvar(3)
		mesh_phys%sos(ix) = sos

		enddo      


	End Subroutine update_phys_Xihua_dudedE

	!> [OK] ener - ein 
	Subroutine update_phys_bkp(ini_cmt,mesh_phys)
		implicit none

		type(inicmt)   :: ini_cmt 
		type(meshphys) :: mesh_phys

		integer :: ix
		real*8  :: den,vel,ene,cvar(3),pvar(3),sos,ein,gamma,vv

		Do ix = 1 , ini_cmt%nx
		!> n ʱ�̵�ֵ
		den = mesh_phys%con_phys(ix,1)
		vel = mesh_phys%con_phys(ix,2) / mesh_phys%con_phys(ix,1)
		ene = mesh_phys%con_phys(ix,3) / mesh_phys%con_phys(ix,1)
		ein = ene - vel*vel/2.d0

		!> n + 1 ʱ�̵�ֵ
		ein = ein - (mesh_phys%pstar(ix+1)*mesh_phys%vel_p(ix+1) - &
			mesh_phys%pstar(ix)  *mesh_phys%vel_p(ix))*ini_cmt%dt / mesh_phys%mass(ix) &
			+ (mesh_phys%pstar(ix+1)*vel - mesh_phys%pstar(ix)*vel)*ini_cmt%dt / mesh_phys%mass(ix) &
			- 0.5d0 * (ini_cmt%dt / mesh_phys%mass(ix))**2 * (mesh_phys%pstar(ix+1)- mesh_phys%pstar(ix))**2
		mesh_phys%ein(ix) = ein

		mesh_phys%con_phys(ix,1) = mesh_phys%mass(ix) / mesh_phys%vol_size_new(ix)
		vel = vel - (mesh_phys%pstar(ix+1)-mesh_phys%pstar(ix))*ini_cmt%dt/mesh_phys%mass(ix)
		ene = ene - (mesh_phys%pstar(ix+1)*mesh_phys%vel_p(ix+1) - &
			mesh_phys%pstar(ix)  *mesh_phys%vel_p(ix))*ini_cmt%dt / mesh_phys%mass(ix)

		vv = (ene - ein) * 2.d0
		if(vv < 1.d-16) then
			vv = 0.d0
		else
			vv = sqrt(vv)
		endif

		if( vel > 0.d0) then
			vel = vv
		else
			vel = -vv
		endif   

		mesh_phys%con_phys(ix,2) = mesh_phys%con_phys(ix,1) * vel
		mesh_phys%con_phys(ix,3) = mesh_phys%con_phys(ix,1) * Ene

		cvar = mesh_phys%con_phys(ix,1:3)

		Call ieos_to_gamma(mesh_phys%eos(ix),gamma)
		Call CVarToPVar(3, CVar,PVar,gamma,mesh_phys%eos(ix))     
		Call PVarToOther(3, PVar,sos,ein,gamma,mesh_phys%eos(ix))

		mesh_phys%vol_size_old(ix) = mesh_phys%vertex(ix+1) - mesh_phys%vertex(ix)
		mesh_phys%pre(ix) = pvar(3)
		mesh_phys%sos(ix) = sos

		enddo      

	End Subroutine update_phys_bkp

end module mdu_fluid_euler_eqn_solver_centered_update
