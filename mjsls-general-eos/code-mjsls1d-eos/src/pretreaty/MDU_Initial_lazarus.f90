!> Last Modified : Sat Feb 25 20:25:34 2023

!> ��ʼ��Module

Module MDU_Initial_lazarus
	use common_inicmt
	use common_meshphys
	use fluid_euler_eqn_basic
Contains


	Subroutine Initial_Mesh_Phys_lazarus(ini_cmt,mesh_phys)
		implicit none

		type(inicmt)   :: ini_cmt 
		type(meshphys) :: mesh_phys

		integer :: ix
		real*8  :: dx,pvar(3),cvar(3),sos,ein,xmid,gamma,vel_r

		ini_cmt%time_stop = 0.8d0 
		!ini_cmt%time_stop = 0.74d0

		dx = 1.d0 / (ini_cmt%nx + 0.d0)

		do ix = 1 , ini_cmt%nx + 1
		mesh_phys%vertex(ix) = (ix - 1) * dx
		enddo
		mesh_phys%vel_p = 0.d0
		mesh_phys%pstar = 0.d0

		Call lazarus_velocity(0.d0,vel_r)

		do ix = 1 , ini_cmt%nx
		xmid = (mesh_phys%vertex(ix) + mesh_phys%vertex(ix+1)) / 2.d0
		PVar(1) = 1.d0      !> rho
		PVar(2) = 0.d0 !vel_r      !> u
		PVar(3) = (5.d0/3.d0-1.d0)*1.d0*1.d-5      !> p
		mesh_phys%eos(ix) = 513 

		Call ieos_to_gamma(mesh_phys%eos(ix),gamma)

		Call PVarToCVar(3, PVar,CVar,gamma,mesh_phys%eos(ix))

		Call PVarToOther(3, PVar,sos,ein,gamma,mesh_phys%eos(ix))

		mesh_phys%con_phys(ix,1:3) = CVar
		mesh_phys%pre(ix) = pvar(3)
		mesh_phys%sos(ix) = sos
		mesh_phys%ein(ix) = 1.d-5

		enddo

		!> 0:Wall, -1:Given Velocity, -5:Given Pressure, -2:InFlow OutFlow
		mesh_phys%bnd_con(1) = 0
		mesh_phys%bnd_con(2) = 0
		mesh_phys%bnd_vel_pre = 0.d0    


	End Subroutine Initial_Mesh_Phys_lazarus

	subroutine lazarus_velocity(t,vel_r)
		implicit none
		real*8,intent(in)  :: t
		real*8,intent(out) :: vel_r
		real*8  :: f

		f = 1.d0 - 0.185d0 * t - 0.28d0 * t**3

		vel_r = -0.6883545d0 * f / (1.d0 - f*t)**(1.d0-0.6883545d0)

	end subroutine lazarus_velocity


End Module MDU_Initial_lazarus
