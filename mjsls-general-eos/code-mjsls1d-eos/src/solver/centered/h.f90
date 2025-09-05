!> Last Modified : 日 10/ 8 17:36:58 2023

module mdu_fluid_euler_eqn_solver_centered_node_solver
	use common_inicmt
	use common_meshphys
	use fluid_euler_eqn_basic
	use MDU_Initial
	use dt
	use output
	use MDU_exact_riemann_slover
	implicit none

contains

!> node_solver 总控
subroutine fluid_euler_eqn_solver_centered_node_solver(ini_cmt,mesh_phys)
	implicit none

	type(inicmt)   :: ini_cmt 
	type(meshphys) :: mesh_phys  

	select case (ini_cmt%Node_Solver)

	case ('Acoustic')

		call cal_centered_node_solver_acoustic(ini_cmt,mesh_phys)

	case ('MJCLS')

		call cal_centered_node_solver_MJCLS(ini_cmt,mesh_phys)

	case ('MJCLS2P2H')

		call cal_centered_node_solver_MJCLS2P2H(ini_cmt,mesh_phys)

	case default

		print*, 'ini_cmt%Node_Solver = ', ini_cmt%Node_Solver
		stop 'ini_cmt%Node_Solver is error, xihua35'

	end select

	call cal_boundary_node_solver(ini_cmt,mesh_phys)

end subroutine fluid_euler_eqn_solver_centered_node_solver


Subroutine cal_centered_node_solver_MJCLS2P2H(ini_cmt,mesh_phys)
	implicit none

	type(inicmt)   :: ini_cmt 
	type(meshphys) :: mesh_phys   

	integer :: ix,iy
	real*8  :: rho1,c1,a1,u1,p1,v1,pp1,aa1
	real*8  :: rho2,c2,a2,u2,p2,v2,pp2,aa2,dx,dx1,dt,dt1,dt2
	real*8  :: beta1,beta2,us,ubar,pstar,tmpa,tmpb,ps1,ps2,gamma
	real*8  :: ustar1,ustar2,aaa,bbb,ccc,dtau1,dtau2,delta1,m1,m2
	real*8  :: ubar1,pstar1

	dt = ini_cmt%dt

	!> calculate inner node
	do ix = 2 , ini_cmt%nx

		call ieos_to_gamma(mesh_phys%eos(ix),gamma)
		rho1 = mesh_phys%con_phys(ix-1,1)
		c1   = mesh_phys%sos(ix-1)
		a1   = mesh_phys%con_phys(ix-1,1) * mesh_phys%sos(ix-1)
		u1   = mesh_phys%con_phys(ix-1,2)/mesh_phys%con_phys(ix-1,1)
		p1   = mesh_phys%pre(ix-1)
		v1   = mesh_phys%vertex(ix) - mesh_phys%vertex(ix-1)
		m1   = mesh_phys%mass(ix-1)
		dt1  = v1/2.d0/mesh_phys%sos(ix-1)

		rho2 = mesh_phys%con_phys(ix,1)
		c2   = mesh_phys%sos(ix)
		a2   = mesh_phys%con_phys(ix,1) * mesh_phys%sos(ix)
		u2   = mesh_phys%con_phys(ix,2)/mesh_phys%con_phys(ix,1)
		p2   = mesh_phys%pre(ix)
		v2   = mesh_phys%vertex(ix+1) - mesh_phys%vertex(ix)
		m2   = mesh_phys%mass(ix)
		dt2 = v2/2.d0/mesh_phys%sos(ix)

		pstar1 = (p1/a1+p2/a2+(u1-u2)) / (1.d0/a1 + 1.d0/a2)

		if ( u1 > u2  ) then 		
			ubar  = (a1*u1+a2*u2)/(a1+a2) + (p1-p2)/(a1+a2)
			pstar = (p1/a1+p2/a2+(u1-u2)) / (1.d0/a1 + 1.d0/a2)   
		else 
			ubar  = (p1*u1+p2*u2)/(p1+p2)
			pstar = (p1+p2)/2.d0 - gamma * (p1+p2)/2.d0 * dt * (u2-u1) / (v1+v2) 
		endif 
		mesh_phys%vel_p(ix) = ubar  + (p1-p2)/(m1+m2)*dt	
		mesh_phys%pstar(ix) = pstar + (gamma+1.d0)/2.d0*(m1+m2)/(v1+v2)*((p1-p2)/(m1+m2)*dt)**2

	enddo   

End Subroutine cal_centered_node_solver_MJCLS2P2H

Subroutine cal_centered_node_solver_MJCLS(ini_cmt,mesh_phys)
	implicit none

	type(inicmt)   :: ini_cmt 
	type(meshphys) :: mesh_phys   

	integer :: ix,iy
	real*8  :: rho1,c1,a1,u1,p1,v1,pp1,aa1
	real*8  :: rho2,c2,a2,u2,p2,v2,pp2,aa2
	real*8  :: beta1,beta2,us,ubar,pstar,tmpa,tmpb,ps1,ps2,gamma
	real*8  :: ustar1,ustar2,aaa,bbb,ccc,dtau1,dtau2,delta1

	!> calculate inner node
	do ix = 2 , ini_cmt%nx

	call ieos_to_gamma(mesh_phys%eos(ix),gamma)
	rho1 = mesh_phys%con_phys(ix-1,1)
	c1   = mesh_phys%sos(ix-1)
	a1   = mesh_phys%con_phys(ix-1,1) * mesh_phys%sos(ix-1)
	u1   = mesh_phys%con_phys(ix-1,2)/mesh_phys%con_phys(ix-1,1)
	p1   = mesh_phys%pre(ix-1)
	v1   = mesh_phys%vertex(ix) - mesh_phys%vertex(ix-1)

	rho2 = mesh_phys%con_phys(ix,1)
	c2   = mesh_phys%sos(ix)
	a2   = mesh_phys%con_phys(ix,1) * mesh_phys%sos(ix)
	u2   = mesh_phys%con_phys(ix,2)/mesh_phys%con_phys(ix,1)
	p2   = mesh_phys%pre(ix)
	v2   = mesh_phys%vertex(ix+1) - mesh_phys%vertex(ix)

	!> ax^2 + bx + c = 0
	aaa = 0.5d0 * (gamma+1.d0) * (rho1-rho2)
	bbb = - (rho1*u1-rho2*u2) * (gamma+1.d0) - (rho1*c1+rho2*c2)
	ccc = 0.5d0 * (gamma+1.d0) * (rho1*u1**2-rho2*u2**2) + &
		rho1*c1*u1 + rho2*c2*u2 + p1 - p2 

	!> Only one solution
	if ( abs(aaa) < 1.d-6 ) then 
		! ubar  = -ccc / bbb
		! pstar = p1 - rho1*c1*(ubar-u1) + 0.5d0 * rho1* (gamma+1.d0)*(ubar-u1)**2         
		! mesh_phys%vel_p(ix) = ubar
		! mesh_phys%pstar(ix) = pstar  
		
				!> acoustic solver
		ubar  = (a1*u1+a2*u2)/(a1+a2) + (p1-p2)/(a1+a2)
		pstar = (p1/a1+p2/a2+(u1-u2)) / (1.d0/a1 + 1.d0/a2)         
		mesh_phys%vel_p(ix) = ubar
		mesh_phys%pstar(ix) = pstar  			

		cycle 
	endif 

	!> Delta = b^2 - 4ac
	delta1 = bbb*bbb - 4.d0 * aaa * ccc
	if ( delta1 > 0.d0 ) then 
		!> Two solutions
		ustar1 = ( -bbb + sqrt(delta1) ) / 2.d0 / aaa 
		ustar2 = ( -bbb - sqrt(delta1) ) / 2.d0 / aaa 
		ubar  = (a1*u1+a2*u2)/(a1+a2) + (p1-p2)/(a1+a2)

		!> get min distance from acoustic solver
		if ( abs(ubar - ustar1) <= abs(ubar - ustar2) ) then 
			!> dS > 0
			if( abs(rho1*c1*(ustar1-u1))**2>= 0.5d0 * rho1* abs((gamma+1.d0)*(ustar1-u1)**3) .and. &
					abs(rho2*c2*(ustar1-u2))**2>= 0.5d0 * rho2* abs((gamma+1.d0)*(ustar1-u2)**3) )  then 
				ps1 = p1 - rho1*c1*(ustar1-u1) + 0.5d0 * rho1* (gamma+1.d0)*(ustar1-u1)**2
				ps2 = p2 + rho2*c2*(ustar1-u2) + 0.5d0 * rho2* (gamma+1.d0)*(ustar1-u2)**2
				mesh_phys%vel_p(ix) = ustar1
				mesh_phys%pstar(ix) = (ps1 + ps2) / 2.d0 
			else 
				!> acoustic solver
				ubar  = (a1*u1+a2*u2)/(a1+a2) + (p1-p2)/(a1+a2)
				pstar = (p1/a1+p2/a2+(u1-u2)) / (1.d0/a1 + 1.d0/a2)         
				mesh_phys%vel_p(ix) = ubar
				mesh_phys%pstar(ix) = pstar   	
			endif 
		else 
			!> dS > 0
			if( abs(rho1*c1*(ustar2-u1))**2>= 0.5d0 * rho1* abs((gamma+1.d0)*(ustar2-u1)**3) .and. &
					abs(rho2*c2*(ustar2-u2))**2 >= 0.5d0 * rho2* abs((gamma+1.d0)*(ustar2-u2)**3) ) then 
				mesh_phys%vel_p(ix) = ustar2
				ps1 = p1 - rho1*c1*(ustar2-u1) + 0.5d0 * rho1* (gamma+1.d0)*(ustar2-u1)**2
				ps2 = p2 + rho2*c2*(ustar2-u2) + 0.5d0 * rho2* (gamma+1.d0)*(ustar2-u2)**2
				mesh_phys%pstar(ix) = (ps1 + ps2) / 2.d0 
			else 
				!> acoustic solver
				ubar  = (a1*u1+a2*u2)/(a1+a2) + (p1-p2)/(a1+a2)
				pstar = (p1/a1+p2/a2+(u1-u2)) / (1.d0/a1 + 1.d0/a2)         
				mesh_phys%vel_p(ix) = ubar
				mesh_phys%pstar(ix) = pstar 
			endif 

		endif
	else 
		!> acoustic solver
		ubar  = (a1*u1+a2*u2)/(a1+a2) + (p1-p2)/(a1+a2)
		pstar = (p1/a1+p2/a2+(u1-u2)) / (1.d0/a1 + 1.d0/a2)         
		mesh_phys%vel_p(ix) = ubar
		mesh_phys%pstar(ix) = pstar        
	endif 
	enddo   

End Subroutine cal_centered_node_solver_MJCLS


Subroutine cal_centered_node_solver_acoustic(ini_cmt,mesh_phys)
	implicit none

	type(inicmt)   :: ini_cmt 
	type(meshphys) :: mesh_phys   

	integer :: ix,iy
	real*8  :: rho1,c1,a1,u1,p1,v1,pp1,aa1
	real*8  :: rho2,c2,a2,u2,p2,v2,pp2,aa2
	real*8  :: beta1,beta2,us,ubar,pstar,tmpa,tmpb,ps1,ps2,gamma

	!> 计算内点
	do ix = 2 , ini_cmt%nx

	call ieos_to_gamma(mesh_phys%eos(ix),gamma)
	rho1 = mesh_phys%con_phys(ix-1,1)
	c1   = mesh_phys%sos(ix-1)
	a1   = mesh_phys%con_phys(ix-1,1) * mesh_phys%sos(ix-1)
	u1   = mesh_phys%con_phys(ix-1,2)/mesh_phys%con_phys(ix-1,1)
	p1   = mesh_phys%pre(ix-1)
	v1   = mesh_phys%vertex(ix) - mesh_phys%vertex(ix-1)

	rho2 = mesh_phys%con_phys(ix,1)
	c2   = mesh_phys%sos(ix)
	a2   = mesh_phys%con_phys(ix,1) * mesh_phys%sos(ix)
	u2   = mesh_phys%con_phys(ix,2)/mesh_phys%con_phys(ix,1)
	p2   = mesh_phys%pre(ix)
	v2   = mesh_phys%vertex(ix+1) - mesh_phys%vertex(ix)

	!> 弱波近似
	ubar  = (a1*u1+a2*u2)/(a1+a2) + (p1-p2)/(a1+a2)
	pstar = (p1/a1+p2/a2+(u1-u2)) / (1.d0/a1 + 1.d0/a2)

	mesh_phys%vel_p(ix) = ubar
	mesh_phys%pstar(ix) = pstar 

	enddo   

End Subroutine cal_centered_node_solver_acoustic


Subroutine cal_boundary_node_solver(ini_cmt,mesh_phys)
	implicit none

	type(inicmt)   :: ini_cmt 
	type(meshphys) :: mesh_phys   

	real*8    :: vel_r,vel1

	!> 0:Wall, -1:Given Velocity, -5:Given Pressure, -2:InFlow OutFlow
	Select Case (mesh_phys%bnd_con(1))
	Case (0)
		mesh_phys%vel_p(1) = 0.d0
		mesh_phys%pstar(1) = mesh_phys%pre(1)
	Case (-1)
		mesh_phys%vel_p(1) = mesh_phys%bnd_vel_pre(1) 
		mesh_phys%pstar(1) = mesh_phys%pre(1)
	Case (-5)
		stop 'Not Ready'
	Case (-2)
		mesh_phys%vel_p(1) = mesh_phys%con_phys(1,2)/mesh_phys%con_phys(1,1) 
		mesh_phys%pstar(1) = mesh_phys%pre(1)
	Case default
		stop '204,xihuafd7a'
	End Select 

	Select Case (mesh_phys%bnd_con(2))
	Case (0)
		mesh_phys%vel_p(ini_cmt%nx+1) = 0.d0
		mesh_phys%pstar(ini_cmt%nx+1) = mesh_phys%pre(ini_cmt%nx)
	Case (-1)
		mesh_phys%vel_p(ini_cmt%nx+1) = mesh_phys%bnd_vel_pre(2) 
		mesh_phys%pstar(ini_cmt%nx+1) = mesh_phys%pre(ini_cmt%nx)
	Case (-5)
		stop 'Not Ready'
	Case (-2)
		mesh_phys%vel_p(ini_cmt%nx+1) = mesh_phys%con_phys(ini_cmt%nx,2) &
			/mesh_phys%con_phys(ini_cmt%nx,1) 
		mesh_phys%pstar(ini_cmt%nx+1) = mesh_phys%pre(ini_cmt%nx)
	Case default
		stop '220,xihuafeed7a'
	End Select 

end subroutine cal_boundary_node_solver   

end module mdu_fluid_euler_eqn_solver_centered_node_solver
