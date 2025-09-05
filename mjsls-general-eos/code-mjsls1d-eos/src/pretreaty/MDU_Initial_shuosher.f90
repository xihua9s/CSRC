!> Last Modified : Sat Feb 25 22:27:50 2023

Module MDU_Initial_shuosher
	use common_inicmt
	use common_meshphys
	use fluid_euler_eqn_basic
Contains


	Subroutine Initial_Mesh_Phys_shuosher(ini_cmt,mesh_phys)
		implicit none

		type(inicmt)   :: ini_cmt 
		type(meshphys) :: mesh_phys

		integer :: ix
		real*8  :: dx,pvar(3),cvar(3),sos,ein,xmid,gamma

		ini_cmt%time_stop = 1.8d0
		dx = 10.d0 / (ini_cmt%nx + 0.d0)

		do ix = 1 , ini_cmt%nx + 1
		mesh_phys%vertex(ix) = (ix - 1) * dx
		if ( mesh_phys%vertex(ix) < 1.d0) then
			mesh_phys%vel_p(ix) = 2.629369d0
		else 
			mesh_phys%vel_p(ix) = 0.d0
		endif 
		enddo
		mesh_phys%pstar = 0.d0

		do ix = 1 , ini_cmt%nx
		xmid = (mesh_phys%vertex(ix) + mesh_phys%vertex(ix+1)) / 2.d0
		if(xmid < 1.d0) then
			PVar(1) = 3.857143d0      !> rho
			PVar(2) = 2.629369d0      !> u
			PVar(3) = 10.333333d0      !> p
			mesh_phys%eos(ix) = 715
		else
			PVar(1) = 1.d0 + 0.2d0*sin(5.d0*xmid)   !> rho
			PVar(2) = 0.d0      !> u
			PVar(3) = 1.d0     !> p
			mesh_phys%eos(ix) = 715
		endif     

		Call ieos_to_gamma(mesh_phys%eos(ix),gamma)

		Call PVarToCVar(3, PVar,CVar,gamma,mesh_phys%eos(ix))

		Call PVarToOther(3, PVar,sos,ein,gamma,mesh_phys%eos(ix))

		mesh_phys%con_phys(ix,1:3) = CVar
		mesh_phys%pre(ix) = pvar(3)
		mesh_phys%sos(ix) = sos
		mesh_phys%ein(ix) = ein

		enddo

		!> 0:Wall, -1:Given Velocity, -5:Given Pressure, -2:InFlow OutFlow
		mesh_phys%bnd_con(1) = -1
		mesh_phys%bnd_con(2) = 0
		mesh_phys%bnd_vel_pre(1) = 2.629369d0 
		mesh_phys%bnd_vel_pre(2) = 0.d0 

	End Subroutine Initial_Mesh_Phys_shuosher


End Module MDU_Initial_shuosher
