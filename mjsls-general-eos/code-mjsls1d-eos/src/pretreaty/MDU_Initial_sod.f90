!> Last Modified : Sat Feb 25 22:33:31 2023

Module MDU_Initial_sod
	use common_inicmt
	use common_meshphys
	use fluid_euler_eqn_basic
Contains


	Subroutine Initial_Mesh_Phys_sod(ini_cmt,mesh_phys)
		implicit none

		type(inicmt)   :: ini_cmt 
		type(meshphys) :: mesh_phys

		integer :: ix
		real*8  :: dx,pvar(3),cvar(3),sos,ein,xmid,gamma

		ini_cmt%time_stop = 0.2d0
		dx = 1.d0 / (ini_cmt%nx + 0.d0)

		do ix = 1 , ini_cmt%nx + 1
		mesh_phys%vertex(ix) = (ix - 1) * dx
		enddo
		mesh_phys%vel_p = 0.d0
		mesh_phys%pstar = 0.d0

		do ix = 1 , ini_cmt%nx
		xmid = (mesh_phys%vertex(ix) + mesh_phys%vertex(ix+1)) / 2.d0
		if(xmid < 0.5d0) then
			PVar(1) = 1.d0      !> rho
			PVar(2) = 0.d0      !> u
			PVar(3) = 1.d0      !> p
			mesh_phys%eos(ix) = 715
		else
			PVar(1) = 0.125d0   !> rho
			PVar(2) = 0.d0      !> u
			PVar(3) = 0.1d0     !> p
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
		mesh_phys%bnd_con(1) = 0
		mesh_phys%bnd_con(2) = 0
		mesh_phys%bnd_vel_pre = 0.d0 

	End Subroutine Initial_Mesh_Phys_sod


End Module MDU_Initial_sod
