
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

		Case('MJCLS2P2H')

			call update_phys_details(ini_cmt,mesh_phys)      

		Case default

			print*, 'ini_cmt%node_solver is', ini_cmt%node_solver
			stop 'ini_cmt%node_solver is Error'

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

end module mdu_fluid_euler_eqn_solver_centered_update
