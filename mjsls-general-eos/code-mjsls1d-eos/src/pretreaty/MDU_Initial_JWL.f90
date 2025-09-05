!> Last Modified : Sat Feb 25 22:33:31 2023

Module MDU_Initial_JWL
	use common_inicmt
	use common_meshphys
	use fluid_euler_eqn_basic
Contains


	Subroutine Initial_Mesh_Phys_CC_33(ini_cmt,mesh_phys)
		implicit none

		type(inicmt)   :: ini_cmt 
		type(meshphys) :: mesh_phys

		integer :: ix
		real*8  :: dx,pvar(3),cvar(3),sos,ein,xmid,gamma

		ini_cmt%time_stop = 40.d-6
		dx = 1.d0 / (ini_cmt%nx + 0.d0)

		do ix = 1 , ini_cmt%nx + 1
			mesh_phys%vertex(ix) = (ix - 1) * dx
		enddo
		mesh_phys%vel_p = 1000.d0
		mesh_phys%pstar = 0.d0

		do ix = 1 , ini_cmt%nx
		xmid = (mesh_phys%vertex(ix) + mesh_phys%vertex(ix+1)) / 2.d0
		if(xmid < 0.5d0) then
			PVar(1) = 1134.d0      !> rho
			PVar(2) = 1000.d0      !> u
			PVar(3) = 20.d9      !> p
			mesh_phys%eos(ix) = 33
		else
			PVar(1) = 500.d0   !> rho
			PVar(2) = 1000.d0      !> u
			PVar(3) = 20.d9     !> p
			mesh_phys%eos(ix) = 33
		endif     


		Call PVarToCVar_CC(3, PVar,CVar,mesh_phys%eos(ix))

		Call PVarToOther_CC(3, PVar,sos,ein,mesh_phys%eos(ix))

		mesh_phys%con_phys(ix,1:3) = CVar
		mesh_phys%pre(ix) = pvar(3)
		mesh_phys%sos(ix) = sos
		mesh_phys%ein(ix) = ein

		enddo

		!> 0:Wall, -1:Given Velocity, -5:Given Pressure, -2:InFlow OutFlow
		mesh_phys%bnd_con(1) = -2
		mesh_phys%bnd_con(2) = -2
		mesh_phys%bnd_vel_pre = 0.d0 

	End Subroutine Initial_Mesh_Phys_CC_33


	Subroutine Initial_Mesh_Phys_JWL_LX17(ini_cmt,mesh_phys)
		implicit none

		type(inicmt)   :: ini_cmt 
		type(meshphys) :: mesh_phys

		integer :: ix
		real*8  :: dx,pvar(3),cvar(3),sos,ein,xmid,gamma

		ini_cmt%time_stop = 2.d-5
		dx = 1.d0 / (ini_cmt%nx + 0.d0)

		do ix = 1 , ini_cmt%nx + 1
			mesh_phys%vertex(ix) = (ix - 1) * dx
		enddo
		mesh_phys%vel_p = 0.d0
		mesh_phys%pstar = 0.d0

		do ix = 1 , ini_cmt%nx
		xmid = (mesh_phys%vertex(ix) + mesh_phys%vertex(ix+1)) / 2.d0
		if(xmid < 0.5d0) then
			PVar(1) = 952.5d0      !> rho
			PVar(2) = 0.d0      !> u
			PVar(3) = 100.d9      !> p
			mesh_phys%eos(ix) = 17
		else
			PVar(1) = 3810.d0   !> rho
			PVar(2) = 0.d0      !> u
			PVar(3) = 200.d9     !> p
			mesh_phys%eos(ix) = 17
		endif     


		Call PVarToCVar_JWL(3, PVar,CVar,mesh_phys%eos(ix))

		Call PVarToOther_JWL(3, PVar,sos,ein,mesh_phys%eos(ix))

		mesh_phys%con_phys(ix,1:3) = CVar
		mesh_phys%pre(ix) = pvar(3)
		mesh_phys%sos(ix) = sos
		mesh_phys%ein(ix) = ein

		enddo

		!> 0:Wall, -1:Given Velocity, -5:Given Pressure, -2:InFlow OutFlow
		mesh_phys%bnd_con(1) = 0
		mesh_phys%bnd_con(2) = 0
		mesh_phys%bnd_vel_pre = 0.d0 

	End Subroutine Initial_Mesh_Phys_JWL_LX17


	Subroutine Initial_Mesh_Phys_JWL_LX17_Blast(ini_cmt,mesh_phys)
		implicit none

		type(inicmt)   :: ini_cmt 
		type(meshphys) :: mesh_phys

		integer :: ix
		real*8  :: dx,pvar(3),cvar(3),sos,ein,xmid,gamma

		ini_cmt%time_stop = 0.038d0
		dx = 1.d0 / (ini_cmt%nx + 0.d0)

		do ix = 1 , ini_cmt%nx + 1
			mesh_phys%vertex(ix) = (ix - 1) * dx
		enddo
		mesh_phys%vel_p = 0.d0
		mesh_phys%pstar = 0.d0

		do ix = 1 , ini_cmt%nx
			xmid = (mesh_phys%vertex(ix) + mesh_phys%vertex(ix+1)) / 2.d0
			if(xmid < 0.1d0) then
				PVar(1) = 1.d0      !> rho
				PVar(2) = 0.d0      !> u
				PVar(3) = 1d3      !> p
				mesh_phys%eos(ix) = 95022
			elseif(xmid < 0.9d0) then 
				PVar(1) = 1.d0   !> rho
				PVar(2) = 0.d0      !> u
				PVar(3) = 1.d-2     !> p
				mesh_phys%eos(ix) = 95022
			else
				PVar(1) = 1d0    !> rho
				PVar(2) = 0.d0      !> u
				PVar(3) = 1d2     !> p
				mesh_phys%eos(ix) = 95022
			endif     

			Call PVarToCVar_JWL(3, PVar,CVar,mesh_phys%eos(ix))

			Call PVarToOther_JWL(3, PVar,sos,ein,mesh_phys%eos(ix))

			mesh_phys%con_phys(ix,1:3) = CVar
			mesh_phys%pre(ix) = pvar(3)
			mesh_phys%sos(ix) = sos
			mesh_phys%ein(ix) = ein

		enddo

		!> 0:Wall, -1:Given Velocity, -5:Given Pressure, -2:InFlow OutFlow
		mesh_phys%bnd_con(1) = 0
		mesh_phys%bnd_con(2) = 0
		mesh_phys%bnd_vel_pre = 0.d0 

	End Subroutine Initial_Mesh_Phys_JWL_LX17_Blast


	Subroutine Initial_Mesh_Phys_JWL_PBX9502(ini_cmt,mesh_phys)
		implicit none

		type(inicmt)   :: ini_cmt 
		type(meshphys) :: mesh_phys

		integer :: ix
		real*8  :: dx,pvar(3),cvar(3),sos,ein,xmid,gamma

		ini_cmt%time_stop = 160.d-6
		dx = 1.d0 / (ini_cmt%nx + 0.d0)

		do ix = 1 , ini_cmt%nx + 1
			mesh_phys%vertex(ix) = (ix - 1) * dx
		enddo
		mesh_phys%vel_p = 0.d0
		mesh_phys%pstar = 0.d0

		do ix = 1 , ini_cmt%nx
			xmid = (mesh_phys%vertex(ix) + mesh_phys%vertex(ix+1)) / 2.d0
			if(xmid < 0.1d0) then
				PVar(1) = 2025.69d0      !> rho
				PVar(2) = 291.74d0      !> u
				PVar(3) = 12.5d9      !> p
				mesh_phys%eos(ix) = 9502
				mesh_phys%vel_p(ix) = 291.74d0
			elseif(xmid < 0.5d0) then 
				PVar(1) = 1895.d0   !> rho
				PVar(2) = 0.d0      !> u
				PVar(3) = 10.d9     !> p
				mesh_phys%eos(ix) = 9502
			else
				PVar(1) = 1895.d0 * 2.d0    !> rho
				PVar(2) = 0.d0      !> u
				PVar(3) = 10.d9     !> p
				mesh_phys%eos(ix) = 9502
			endif     

			Call PVarToCVar_JWL(3, PVar,CVar,mesh_phys%eos(ix))

			Call PVarToOther_JWL(3, PVar,sos,ein,mesh_phys%eos(ix))

			mesh_phys%con_phys(ix,1:3) = CVar
			mesh_phys%pre(ix) = pvar(3)
			mesh_phys%sos(ix) = sos
			mesh_phys%ein(ix) = ein

		enddo

		!> 0:Wall, -1:Given Velocity, -5:Given Pressure, -2:InFlow OutFlow
		mesh_phys%bnd_con(1) = -2
		mesh_phys%bnd_con(2) = 0
		mesh_phys%bnd_vel_pre = 0.d0 

	End Subroutine Initial_Mesh_Phys_JWL_PBX9502


	Subroutine Initial_Mesh_Phys_JWL_PBX9502_Sod(ini_cmt,mesh_phys)
		implicit none

		type(inicmt)   :: ini_cmt 
		type(meshphys) :: mesh_phys

		integer :: ix
		real*8  :: dx,pvar(3),cvar(3),sos,ein,xmid,gamma

		ini_cmt%time_stop = 1d-5
		dx = 1.d0 / (ini_cmt%nx + 0.d0)

		do ix = 1 , ini_cmt%nx + 1
			mesh_phys%vertex(ix) = (ix - 1) * dx
		enddo
		mesh_phys%vel_p = 0.d0
		mesh_phys%pstar = 0.d0

		do ix = 1 , ini_cmt%nx
			xmid = (mesh_phys%vertex(ix) + mesh_phys%vertex(ix+1)) / 2.d0
			if(xmid < 0.43d0) then
				PVar(1) = 1895.d0      !> rho
				PVar(2) = 0.d0      !> u
				PVar(3) = 10d11      !> p
				mesh_phys%eos(ix) = 95021
				mesh_phys%vel_p(ix) = 0d0
			else
				PVar(1) = 236.875d0 !> rho
				PVar(2) = 0.d0      !> u
				PVar(3) = 10.d10     !> p
				mesh_phys%eos(ix) = 95021
				mesh_phys%vel_p(ix) = 0d0
			endif     

			Call PVarToCVar_JWL(3, PVar,CVar,mesh_phys%eos(ix))

			Call PVarToOther_JWL(3, PVar,sos,ein,mesh_phys%eos(ix))

			mesh_phys%con_phys(ix,1:3) = CVar
			mesh_phys%pre(ix) = pvar(3)
			mesh_phys%sos(ix) = sos
			mesh_phys%ein(ix) = ein

		enddo

		!> 0:Wall, -1:Given Velocity, -5:Given Pressure, -2:InFlow OutFlow
		mesh_phys%bnd_con(1) = 0
		mesh_phys%bnd_con(2) = 0
		mesh_phys%bnd_vel_pre = 0.d0 

	End Subroutine Initial_Mesh_Phys_JWL_PBX9502_Sod



	Subroutine Initial_Mesh_Phys_vandeW_Sod(ini_cmt,mesh_phys)
		implicit none

		type(inicmt)   :: ini_cmt 
		type(meshphys) :: mesh_phys

		integer :: ix
		real*8  :: dx,pvar(3),cvar(3),sos,ein,xmid,gamma

		ini_cmt%time_stop = 0.25d0
		dx = 1.d0 / (ini_cmt%nx + 0.d0)

		do ix = 1 , ini_cmt%nx + 1
			mesh_phys%vertex(ix) = (ix - 1) * dx
		enddo
		mesh_phys%vel_p = 0.d0
		mesh_phys%pstar = 0.d0

		do ix = 1 , ini_cmt%nx
			xmid = (mesh_phys%vertex(ix) + mesh_phys%vertex(ix+1)) / 2.d0
			if(xmid < 0.4d0) then
				PVar(1) = 1.d0      !> rho
				PVar(2) = 0.d0      !> u
				PVar(3) = 1.d0       !> p
				mesh_phys%eos(ix) = 22423
				mesh_phys%vel_p(ix) = 0d0
			else
				PVar(1) = 0.1d0 !> rho
				PVar(2) = 0.d0      !> u
				PVar(3) = 0.125d0     !> p
				mesh_phys%eos(ix) = 22423
				mesh_phys%vel_p(ix) = 0d0
			endif     

			Call PVarToCVar_vandeW(3, PVar,CVar,mesh_phys%eos(ix))

			Call PVarToOther_vandeW(3, PVar,sos,ein,mesh_phys%eos(ix))

			mesh_phys%con_phys(ix,1:3) = CVar
			mesh_phys%pre(ix) = pvar(3)
			mesh_phys%sos(ix) = sos
			mesh_phys%ein(ix) = ein

		enddo

		!> 0:Wall, -1:Given Velocity, -5:Given Pressure, -2:InFlow OutFlow
		mesh_phys%bnd_con(1) = 0
		mesh_phys%bnd_con(2) = 0
		mesh_phys%bnd_vel_pre = 0.d0 

	End Subroutine Initial_Mesh_Phys_vandeW_Sod




	Subroutine Initial_Mesh_Phys_Wilkins(ini_cmt,mesh_phys)
		implicit none

		type(inicmt)   :: ini_cmt 
		type(meshphys) :: mesh_phys

		integer :: ix
		real*8  :: dx,pvar(3),cvar(3),sos,ein,xmid,gamma

		ini_cmt%time_stop = 5d-6
		dx = 0.05d0 / (ini_cmt%nx + 0.d0)

		do ix = 1 , ini_cmt%nx + 1
			mesh_phys%vertex(ix) = (ix - 1) * dx
		enddo
		mesh_phys%vel_p = 0.d0
		mesh_phys%pstar = 0.d0

		do ix = 1 , ini_cmt%nx
			xmid = (mesh_phys%vertex(ix) + mesh_phys%vertex(ix+1)) / 2.d0
			if(xmid < 0.005d0) then
				PVar(1) = 2785.d0      !> rho
				PVar(2) = 800.d0      !> u
				PVar(3) = 1.d-6      !> p
				mesh_phys%eos(ix) = 14719
				mesh_phys%vel_p(ix) = 800d0
			else
				PVar(1) = 2785.d0      !> rho
				PVar(2) = 0.d0      !> u
				PVar(3) = 1.d-6      !> p
				mesh_phys%eos(ix) = 14719
				mesh_phys%vel_p(ix) = 0.d0
			endif     

			Call PVarToCVar_MieGruSolid(3, PVar,CVar,mesh_phys%eos(ix))

			Call PVarToOther_MieGruSolid(3, PVar,sos,ein,mesh_phys%eos(ix))

			mesh_phys%con_phys(ix,1:3) = CVar
			mesh_phys%pre(ix) = pvar(3)
			mesh_phys%sos(ix) = sos
			mesh_phys%ein(ix) = ein

		enddo

		!> 0:Wall, -1:Given Velocity, -5:Given Pressure, -2:InFlow OutFlow
		mesh_phys%bnd_con(1) = -2
		mesh_phys%bnd_con(2) = 0
		mesh_phys%bnd_vel_pre = 1.d-6

	End Subroutine Initial_Mesh_Phys_Wilkins

End Module MDU_Initial_JWL
