!> Last Modified : 日  4/23 21:46:26 2023

Module MDU_Initial_sedov
	use common_inicmt
	use common_meshphys
	use fluid_euler_eqn_basic
Contains

	Subroutine Initial_Mesh_Phys_sedov_2(ini_cmt,mesh_phys)
		implicit none

		type(inicmt)   :: ini_cmt 
		type(meshphys) :: mesh_phys

		integer :: ix
		real*8  :: dx,pvar(3),cvar(3),sos,ein,xmid,gamma,volume

		ini_cmt%time_stop = 1.d-3
		dx = 4d0 / (ini_cmt%nx + 0.d0)

		do ix = 1 , ini_cmt%nx + 1
		mesh_phys%vertex(ix) = (ix - 1) * dx
		enddo
		mesh_phys%vel_p = 0.d0
		mesh_phys%pstar = 0.d0

		do ix = 1 , ini_cmt%nx
		xmid = (mesh_phys%vertex(ix) + mesh_phys%vertex(ix+1)) / 2.d0
		PVar(1) = 1.d0      !> rho
		PVar(2) = 0.d0      !> u
		PVar(3) = (7.d0/5.d0-1.d0)*1.d0*1.d-10      !> p
		mesh_phys%eos(ix) = 715   

		if(ix == ini_cmt%nx/2 .or. ix == ini_cmt%nx/2+1) then
			select case (ini_cmt%cacysp) 
			case('Cartesian')
				mesh_phys%ein(ix) =  1600000.d0 / dx                   
				mesh_phys%con_phys(ix,3) = mesh_phys%ein(ix)
				pvar(3) = (7.d0/5.d0-1.d0)*mesh_phys%ein(ix)
			case('Cylindrical')
				mesh_phys%ein(ix) =  0.306831354d0 / dx /dx          
				mesh_phys%con_phys(ix,3) = mesh_phys%ein(ix)
				pvar(3) = (7.d0/5.d0-1.d0)*mesh_phys%ein(ix)
			case('Spheral')
				mesh_phys%ein(ix) =  0.244816d0 * (3.14159265359/4.d0) / dx**3          
				mesh_phys%con_phys(ix,3) = mesh_phys%ein(ix)
				pvar(3) = (7.d0/5.d0-1.d0)*mesh_phys%ein(ix)
			case default
				print*, 'ini_cmt%cacysp=', ini_cmt%cacysp
				stop 'ini_cmt%cacysp is error.'
			end select
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


	End Subroutine Initial_Mesh_Phys_sedov_2


	Subroutine Initial_Mesh_Phys_sedov(ini_cmt,mesh_phys)
		implicit none

		type(inicmt)   :: ini_cmt 
		type(meshphys) :: mesh_phys

		integer :: ix
		real*8  :: dx,pvar(3),cvar(3),sos,ein,xmid,gamma,volume

		ini_cmt%time_stop = 1.d0
		dx = 1.2d0 / (ini_cmt%nx + 0.d0)

		do ix = 1 , ini_cmt%nx + 1
		mesh_phys%vertex(ix) = (ix - 1) * dx
		enddo
		mesh_phys%vel_p = 0.d0
		mesh_phys%pstar = 0.d0

		do ix = 1 , ini_cmt%nx
		xmid = (mesh_phys%vertex(ix) + mesh_phys%vertex(ix+1)) / 2.d0
		PVar(1) = 1.d0      !> rho
		PVar(2) = 0.d0      !> u
		PVar(3) = 1.d-6      !> p
		mesh_phys%eos(ix) = 715   

		if(ix == 1) then
			select case (ini_cmt%cacysp) 
			case('Cartesian')
				mesh_phys%ein(ix) =  0.5385952d0 / dx                   
				mesh_phys%con_phys(ix,3) = mesh_phys%ein(ix)
				pvar(3) = (7.d0/5.d0-1.d0)*mesh_phys%ein(ix)
			case('Cylindrical')
				mesh_phys%ein(ix) =  0.306831354d0 / dx /dx          
				mesh_phys%con_phys(ix,3) = mesh_phys%ein(ix)
				pvar(3) = (7.d0/5.d0-1.d0)*mesh_phys%ein(ix)
			case('Spheral')
				mesh_phys%ein(ix) =  0.244816d0 * (3.14159265359/4.d0) / dx**3          
				mesh_phys%con_phys(ix,3) = mesh_phys%ein(ix)
				pvar(3) = (7.d0/5.d0-1.d0)*mesh_phys%ein(ix)
			case default
				print*, 'ini_cmt%cacysp=', ini_cmt%cacysp
				stop 'ini_cmt%cacysp is error.'
			end select
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


	End Subroutine Initial_Mesh_Phys_sedov


End Module MDU_Initial_sedov
