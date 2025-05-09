
module mdu_fluid_euler_eqn_solver_staggered   
use MDU_Initial
use common_inicmt
use common_meshphys
use fluid_euler_eqn_basic
use dt
use output
implicit none

contains

Subroutine fluid_euler_eqn_solver_staggered(ini_cmt,mesh_phys)
	implicit none

	type(inicmt)   :: ini_cmt 
	type(meshphys) :: mesh_phys

	select case (ini_cmt%Node_Solver)

	case ('MJSLS')

		call solver_staggered_MJSLS(ini_cmt,mesh_phys)

	case default

		stop 'Error'

	end select

end Subroutine fluid_euler_eqn_solver_staggered   

Subroutine solver_staggered_MJSLS(ini_cmt,mesh_phys)
implicit none

type(inicmt)   :: ini_cmt 
type(meshphys) :: mesh_phys

Do While (ini_cmt%current_time < ini_cmt%Time_stop)

	call cal_dt(ini_cmt,mesh_phys)

	Call solver_staggered_MJSLS_details( &
		ini_cmt%nx,      &
		ini_cmt%dt,      &
		mesh_phys%con_phys(:,1), &
		mesh_phys%pre, &
		mesh_phys%sos, &
		mesh_phys%ein, &
		mesh_phys%mass, &
		mesh_phys%eos, &
		mesh_phys%vol, &
		mesh_phys%vertex, &
		mesh_phys%vel_p, &
		mesh_phys%bnd_con, &
		mesh_phys%bnd_vel_pre &
		)

	if( mod(ini_cmt%it,ini_cmt%It_out) ==0 ) then
		Call output_mesh_phys(ini_cmt,mesh_phys)
	endif

Enddo   

	call output_mesh_phys_final(ini_cmt,mesh_phys)

end Subroutine solver_staggered_MJSLS


Subroutine solver_staggered_MJSLS_details(nx,dt,den,pre,sos,ein,mass,eos,&
		vol,x,vel,bnd_con,bnd_vel_pre)
	implicit none

	integer,intent(in)   :: nx,bnd_con(2),eos(nx)
	real*8,intent(in)    :: mass(nx),dt,bnd_vel_pre(4)
	real*8,intent(inout) :: den(nx),pre(nx),sos(nx),ein(nx),vol(nx)
	real*8,intent(inout) :: x(nx+1),vel(nx+1)

	integer :: ix
	real*8  :: vis1,vis2,gamma,tau(nx)
	real*8  :: q(nx),vel_pre(nx+1),x_pre(nx+1),pre_pre(nx)
	real*8  :: vol_pre(nx),tau_pre(nx),den_pre(nx),ein_pre(nx),sos_pre(nx)
	real*8  :: q_pre(nx),vel_cor(nx+1),x_cor(nx+1),pre_cor(nx)
	real*8  :: vol_cor(nx),tau_cor(nx),den_cor(nx),ein_cor(nx),sos_cor(nx)
	real*8  :: tau1(nx),pre1(nx),tau2(nx),pre2(nx)

	vis1  = 1.d0
	vis2  = 1.d0    
	call ieos_to_gamma(eos(1),gamma)

	do ix = 1 , nx 
	tau(ix)  = 1.d0 / den(ix)
	tau1(ix) = tau(ix) + (vel(ix+1)-vel(ix))/den(ix)/sos(ix)
	if( vel(ix+1) < vel(ix) ) then
		pre1(ix) = pre(ix) - vis1*(den(ix)*sos(ix))**2*(tau1(ix)-tau(ix)) + &
			vis2*(gamma+1.d0)/2.d0*den(ix)**3*sos(ix)**2*(tau1(ix)-tau(ix))**2 
	else 
		pre1(ix) = pre(ix) 
	endif 
	enddo 

	do ix = 2 , nx 
		vel_pre(ix) = vel(ix) - dt*(pre1(ix) - pre1(ix-1)) / (mass(ix)+mass(ix-1))*2.d0
	enddo 

	select case (bnd_con(1))
	case(0)
		vel_pre(1) = 0.d0 
	case(-1)
		vel_pre(1) = bnd_vel_pre(1) 
	case(-5)
		stop 'Not Ready'
	case(-2)
		vel_pre(1) = vel_pre(2)      
	case default 
		stop 'error329'
	end select 
	select case (bnd_con(2))
	case(0)
		vel_pre(nx+1) = 0.d0 
	case(-1)
		vel_pre(nx+1) = bnd_vel_pre(2) 
	case(-5)
		stop 'Not Ready'
	case(-2)
		vel_pre(nx+1) = vel_pre(nx)      
	case default 
		stop 'error,339'
	end select 

	do ix = 1 , nx + 1
	x_pre(ix) = x(ix) + dt*(vel(ix) + vel_pre(ix))/2.d0
	enddo 

	do ix = 1 , nx 
	vol_pre(ix) = x_pre(ix+1) - x_pre(ix)
	tau_pre(ix) = vol_pre(ix) / mass(ix)
	den_pre(ix) = 1.d0 / tau_pre(ix)
	ein_pre(ix) = ein(ix) - (pre1(ix))*(vol_pre(ix)-vol(ix))/mass(ix)
	pre_pre(ix) = (gamma - 1.d0) * den_pre(ix) * ein_pre(ix)
	sos_pre(ix) = sqrt(gamma*pre_pre(ix)/den_pre(ix))
	enddo 

	do ix = 1 , nx 
		tau2(ix) = tau_pre(ix) + (vel_pre(ix+1)-vel_pre(ix))/den_pre(ix)/sos_pre(ix)
		if( vel_pre(ix+1) < vel_pre(ix) ) then
			pre2(ix) = pre_pre(ix) - vis1*(den_pre(ix)*sos_pre(ix))**2*(tau2(ix)-tau_pre(ix)) + &
				vis2*(gamma+1.d0)/2.d0*den_pre(ix)**3*sos_pre(ix)**2*(tau2(ix)-tau_pre(ix))**2                   
		else
			pre2(ix) = pre_pre(ix) 
		endif 
	enddo

	do ix = 2 , nx 
	vel_cor(ix) = vel(ix) - dt*( pre2(ix) - pre2(ix-1) + pre1(ix) - pre1(ix-1) ) &
		/ (mass(ix)+mass(ix-1))*2.d0 / 2.d0
	enddo 

	select case (bnd_con(1))
	case(0)
		vel_cor(1) = 0.d0 
	case(-1)
		vel_cor(1) = bnd_vel_pre(1) 
	case(-5)
		stop 'Not Ready'
	case(-2)
		vel_cor(1) = vel_cor(2)      
	case default 
		stop 'error59'
	end select 
	select case (bnd_con(2))
	case(0)
		vel_cor(nx+1) = 0.d0 
	case(-1)
		vel_cor(nx+1) = bnd_vel_pre(2) 
	case(-5)
		stop 'Not Ready'
	case(-2)
		vel_cor(nx+1) = vel_cor(nx)      
	case default 
		stop 'error989'
	end select 

	do ix = 1 , nx + 1
		x_cor(ix) = x(ix) + dt* (vel_cor(ix) + vel(ix))/2.d0 
	enddo 

	do ix = 1 , nx 
		vol_cor(ix) = x_cor(ix+1) - x_cor(ix)
		tau_cor(ix) = vol_cor(ix) / mass(ix)
		den_cor(ix) = 1.d0 / tau_cor(ix)
		ein_cor(ix) = ein(ix) - (pre2(ix)+pre1(ix))*(vol_cor(ix)-vol(ix))/mass(ix)/2.d0
		pre_cor(ix) = (gamma - 1.d0)*den_cor(ix)*ein_cor(ix)
		sos_cor(ix) = sqrt(gamma*pre_cor(ix)/den_cor(ix))
	enddo 

	vol = vol_cor
	den = den_cor
	ein = ein_cor
	pre = pre_cor 
	sos = sos_cor
	x   = x_cor
	vel = vel_cor 

end Subroutine solver_staggered_MJSLS_details


end module mdu_fluid_euler_eqn_solver_staggered  
