
module mdu_MJSLS_dpdtau   
use MDU_Initial
use common_inicmt
use common_meshphys
use fluid_euler_eqn_basic
use dt
use output
implicit none

contains

!> Xihua Xu, A parameter-free staggered-grid Lagrangian 
!> scheme for two-dimensional compressible flow problems
!> Journal of Computational Physics, 2024

Subroutine solver_staggered_MJSLSdpdtau(ini_cmt,mesh_phys)
	implicit none

	type(inicmt)   :: ini_cmt 
	type(meshphys) :: mesh_phys

	Do While (ini_cmt%current_time < ini_cmt%Time_stop)

		call cal_dt(ini_cmt,mesh_phys)

		Call solver_staggered_MJSLSdpdtau_details( &
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

end Subroutine solver_staggered_MJSLSdpdtau


Subroutine solver_staggered_MJSLSdpdtau_details(nx,dt,den,pre,sos,ein,mass,eos,&
		vol,x,vel,bnd_con,bnd_vel_pre)
	implicit none

	integer,intent(in)   :: nx,bnd_con(2),eos(nx)
	real*8,intent(in)    :: mass(nx),dt,bnd_vel_pre(4)
	real*8,intent(inout) :: den(nx),pre(nx),sos(nx),ein(nx),vol(nx)
	real*8,intent(inout) :: x(nx+1),vel(nx+1)

	integer :: ix
	real*8  :: vis1,vis2,gamma6(6),tau(nx)
	real*8  :: q(nx),vel_pre(nx+1),x_pre(nx+1),pre_pre(nx)
	real*8  :: vol_pre(nx),tau_pre(nx),den_pre(nx),ein_pre(nx),sos_pre(nx)
	real*8  :: q_pre(nx),vel_cor(nx+1),x_cor(nx+1),pre_cor(nx)
	real*8  :: vol_cor(nx),tau_cor(nx),den_cor(nx),ein_cor(nx),sos_cor(nx)
	real*8  :: tau1(nx),pre1(nx),tau2(nx),pre2(nx),dpdtau,dp2dtau2
  
	do ix = 1 , nx 
		tau(ix)  = 1.d0 / den(ix)
		tau1(ix) = tau(ix) + (vel(ix+1)-vel(ix))/den(ix)/sos(ix)
		call ieos_to_dpdtau(eos(ix),tau(ix),pre(ix),dpdtau,dp2dtau2)
		if( vel(ix+1) < vel(ix) ) then
			pre1(ix) = pre(ix) + dpdtau*(tau1(ix)-tau(ix)) + &
								0.5d0*dp2dtau2*(tau1(ix)-tau(ix))**2 
		else 
			pre1(ix) = pre(ix) 
		endif 
	enddo 

	do ix = 2 , nx 
		vel_pre(ix) = vel(ix) - dt*(pre1(ix) - pre1(ix-1)) / (mass(ix)+mass(ix-1))*2.d0
	enddo 

	!> 0:Wall, -1:Given Velocity, -5:Given Pressure, -2:InFlow OutFlow
	select case (bnd_con(1))
	case(0)
		vel_pre(1) = 0.d0 
	case(-1)
		vel_pre(1) = bnd_vel_pre(1) 
	case(-5)
		vel_pre(1) = bnd_vel_pre(1) 
	case(-2)
		vel_pre(1) = vel_pre(2)      
	case default 
		stop 'errorxihua,329'
	end select 
	select case (bnd_con(2))
	case(0)
		vel_pre(nx+1) = 0.d0 
	case(-1)
		vel_pre(nx+1) = bnd_vel_pre(2) 
	case(-5)
		vel_pre(nx+1) = bnd_vel_pre(2) 
		! stop 'Not Ready'
	case(-2)
		vel_pre(nx+1) = vel_pre(nx)      
	case default 
		stop 'errorxihua,329'
	end select 

	do ix = 1 , nx + 1
		x_pre(ix) = x(ix) + dt*(vel(ix) + vel_pre(ix))/2.d0
	enddo 

	do ix = 1 , nx 
		vol_pre(ix) = x_pre(ix+1) - x_pre(ix)
		tau_pre(ix) = vol_pre(ix) / mass(ix)
		den_pre(ix) = 1.d0 / tau_pre(ix)
		ein_pre(ix) = ein(ix) - (pre1(ix))*(vol_pre(ix)-vol(ix))/mass(ix)
		call den_ein_2_pre_sos(eos(ix),den_pre(ix),ein_pre(ix),pre_pre(ix),sos_pre(ix))
	enddo 

	! vol = vol_pre
	! den = den_pre
	! ein = ein_pre
	! pre = pre_pre 
	! sos = sos_pre
	! x   = x_pre
	! vel = vel_pre
	! return

	do ix = 1 , nx 
		tau2(ix) = tau_pre(ix) + (vel_pre(ix+1)-vel_pre(ix))/den_pre(ix)/sos_pre(ix)
		call ieos_to_dpdtau(eos(ix),tau_pre(ix),pre_pre(ix),dpdtau,dp2dtau2)
		if( vel_pre(ix+1) < vel_pre(ix) ) then
			pre2(ix) = pre_pre(ix) + dpdtau*(tau2(ix)-tau_pre(ix)) + &
							0.5d0*dp2dtau2*(tau2(ix)-tau_pre(ix))**2                   
		else
			pre2(ix) = pre_pre(ix) 
		endif 
	enddo

	do ix = 2 , nx 
		vel_cor(ix) = vel(ix) - dt*( pre2(ix) - pre2(ix-1) + pre1(ix) - pre1(ix-1) ) &
		/ (mass(ix)+mass(ix-1))*2.d0 / 2.d0
	enddo 

	!> 0:Wall, -1:Given Velocity, -5:Given Pressure, -2:InFlow OutFlow
	select case (bnd_con(1))
	case(0)
		vel_cor(1) = 0.d0 
	case(-1)
		vel_cor(1) = bnd_vel_pre(1) 
	case(-5)
		vel_cor(1) = bnd_vel_pre(1) 
		! stop 'Not Ready'
	case(-2)
		vel_cor(1) = vel_cor(2)      
	case default 
		stop 'errorxihua,329'
	end select 
	select case (bnd_con(2))
	case(0)
		vel_cor(nx+1) = 0.d0 
	case(-1)
		vel_cor(nx+1) = bnd_vel_pre(2) 
	case(-5)
		! stop 'Not Ready'
		vel_cor(nx+1) = bnd_vel_pre(2) 
	case(-2)
		vel_cor(nx+1) = vel_cor(nx)      
	case default 
		stop 'errorxihua,329'
	end select 

	do ix = 1 , nx + 1
		x_cor(ix) = x(ix) + dt* (vel_cor(ix) + vel(ix))/2.d0 
	enddo 

	do ix = 1 , nx 
		vol_cor(ix) = x_cor(ix+1) - x_cor(ix)
		tau_cor(ix) = vol_cor(ix) / mass(ix)
		den_cor(ix) = 1.d0 / tau_cor(ix)
		ein_cor(ix) = ein(ix) - (pre2(ix)+pre1(ix))*(vol_cor(ix)-vol(ix))/mass(ix)/2.d0
		call den_ein_2_pre_sos(eos(ix),den_cor(ix),ein_cor(ix),pre_cor(ix),sos_cor(ix))
	enddo 

	vol = vol_cor
	den = den_cor
	ein = ein_cor
	pre = pre_cor 
	sos = sos_cor
	x   = x_cor
	vel = vel_cor 

end Subroutine solver_staggered_MJSLSdpdtau_details


!> for paper

Subroutine solver_staggered_MJSLSdpdtauDT(ini_cmt,mesh_phys)
	implicit none

	type(inicmt)   :: ini_cmt 
	type(meshphys) :: mesh_phys

	Do While (ini_cmt%current_time < ini_cmt%Time_stop)

		call cal_dt(ini_cmt,mesh_phys)

		Call solver_staggered_MJSLSdpdtauDT_details( &
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

end Subroutine solver_staggered_MJSLSdpdtauDT


Subroutine solver_staggered_MJSLSdpdtauDT_details(nx,dt,den,pre,sos,ein,mass,eos,&
		vol,x,vel,bnd_con,bnd_vel_pre)
	implicit none

	integer,intent(in)   :: nx,bnd_con(2),eos(nx)
	real*8,intent(in)    :: mass(nx),dt,bnd_vel_pre(4)
	real*8,intent(inout) :: den(nx),pre(nx),sos(nx),ein(nx),vol(nx)
	real*8,intent(inout) :: x(nx+1),vel(nx+1)

	integer :: ix
	real*8  :: vol1,vol2,gamma6(6),tau(nx),P1,P2
	real*8  :: q(nx),vel_pre(nx+1),x_pre(nx+1),pre_pre(nx)
	real*8  :: vol_pre(nx),tau_pre(nx),den_pre(nx),ein_pre(nx),sos_pre(nx)
	real*8  :: q_pre(nx),vel_cor(nx+1),x_cor(nx+1),pre_cor(nx)
	real*8  :: vol_cor(nx),tau_cor(nx),den_cor(nx),ein_cor(nx),sos_cor(nx)
	real*8  :: tau1(nx),pre1(nx),tau2(nx),pre2(nx),dpdtau,dp2dtau2
  
	do ix = 1 , nx 
		tau(ix)  = 1.d0 / den(ix)
		call ieos_to_dpdtau(eos(ix),tau(ix),pre(ix),dpdtau,dp2dtau2)

		if( vel(ix+1) < vel(ix) ) then

			vol1 = (x(ix+1)-x(ix))-(vel(ix+1)-vel(ix))*(x(ix+1)-x(ix))/sos(ix)
			vol2 = (x(ix+1)-x(ix))-(vel(ix+1)-vel(ix))*dt

			P1 = pre(ix) + dpdtau*((vel(ix+1)-vel(ix))/den(ix)/sos(ix)) + &
					 0.5d0*dp2dtau2*((vel(ix+1)-vel(ix))/den(ix)/sos(ix))**2
			P2 = pre(ix) + dpdtau*((vel(ix+1)-vel(ix))*dt/den(ix)/(x(ix+1)-x(ix))) + &
					 0.5d0*dp2dtau2*((vel(ix+1)-vel(ix))*dt/den(ix)/(x(ix+1)-x(ix)))**2

			pre1(ix) = (P1*vol1+P2*vol2)/(vol1+vol2)
		else 
			pre1(ix) = pre(ix) 
		endif 
	enddo 

	do ix = 2 , nx 
		vel_pre(ix) = vel(ix) - dt*(pre1(ix) - pre1(ix-1)) / (mass(ix)+mass(ix-1))*2.d0
	enddo 

	!> 0:Wall, -1:Given Velocity, -5:Given Pressure, -2:InFlow OutFlow
	select case (bnd_con(1))
	case(0)
		vel_pre(1) = 0.d0 
	case(-1)
		vel_pre(1) = bnd_vel_pre(1) 
	case(-5)
		vel_pre(1) = bnd_vel_pre(1) 
	case(-2)
		vel_pre(1) = vel_pre(2)      
	case default 
		stop 'errorxihua,329'
	end select 
	select case (bnd_con(2))
	case(0)
		vel_pre(nx+1) = 0.d0 
	case(-1)
		vel_pre(nx+1) = bnd_vel_pre(2) 
	case(-5)
		vel_pre(nx+1) = bnd_vel_pre(2) 
		! stop 'Not Ready'
	case(-2)
		vel_pre(nx+1) = vel_pre(nx)      
	case default 
		stop 'errorxihua,329'
	end select 

	do ix = 1 , nx + 1
		x_pre(ix) = x(ix) + dt*(vel(ix) + vel_pre(ix))/2.d0
	enddo 

	do ix = 1 , nx 
		vol_pre(ix) = x_pre(ix+1) - x_pre(ix)
		tau_pre(ix) = vol_pre(ix) / mass(ix)
		den_pre(ix) = 1.d0 / tau_pre(ix)
		ein_pre(ix) = ein(ix) - (pre1(ix))*(vol_pre(ix)-vol(ix))/mass(ix)
		call den_ein_2_pre_sos(eos(ix),den_pre(ix),ein_pre(ix),pre_pre(ix),sos_pre(ix))
	enddo 

	! vol = vol_pre
	! den = den_pre
	! ein = ein_pre
	! pre = pre_pre 
	! sos = sos_pre
	! x   = x_pre
	! vel = vel_pre
	! return

	do ix = 1 , nx 
 
		call ieos_to_dpdtau(eos(ix),tau_pre(ix),pre_pre(ix),dpdtau,dp2dtau2)

		if( vel_pre(ix+1) < vel_pre(ix) ) then

			vol1 = (x_pre(ix+1)-x_pre(ix))-(vel_pre(ix+1)-vel_pre(ix))*&
					(x_pre(ix+1)-x_pre(ix))/sos_pre(ix)
			vol2 = (x_pre(ix+1)-x_pre(ix))-(vel_pre(ix+1)-vel_pre(ix))*dt

			P1 = pre_pre(ix) + dpdtau*((vel_pre(ix+1)-vel_pre(ix))/den_pre(ix)/sos_pre(ix)) + &
					 0.5d0*dp2dtau2*((vel_pre(ix+1)-vel_pre(ix))/den_pre(ix)/sos_pre(ix))**2
			P2 = pre_pre(ix) + dpdtau*((vel_pre(ix+1)-vel_pre(ix))*dt/ &
										den_pre(ix)/(x_pre(ix+1)-x_pre(ix))) + &
					 0.5d0*dp2dtau2*((vel_pre(ix+1)-vel_pre(ix))*dt &
					 					/den_pre(ix)/(x_pre(ix+1)-x_pre(ix)))**2

			pre2(ix) = (P1*vol1+P2*vol2)/(vol1+vol2)
		else
			pre2(ix) = pre_pre(ix) 
		endif 
	enddo

	do ix = 2 , nx 
		vel_cor(ix) = vel(ix) - dt*( pre2(ix) - pre2(ix-1) + pre1(ix) - pre1(ix-1) ) &
		/ (mass(ix)+mass(ix-1))*2.d0 / 2.d0
	enddo 

	!> 0:Wall, -1:Given Velocity, -5:Given Pressure, -2:InFlow OutFlow
	select case (bnd_con(1))
	case(0)
		vel_cor(1) = 0.d0 
	case(-1)
		vel_cor(1) = bnd_vel_pre(1) 
	case(-5)
		vel_cor(1) = bnd_vel_pre(1) 
		! stop 'Not Ready'
	case(-2)
		vel_cor(1) = vel_cor(2)      
	case default 
		stop 'errorxihua,329'
	end select 
	select case (bnd_con(2))
	case(0)
		vel_cor(nx+1) = 0.d0 
	case(-1)
		vel_cor(nx+1) = bnd_vel_pre(2) 
	case(-5)
		! stop 'Not Ready'
		vel_cor(nx+1) = bnd_vel_pre(2) 
	case(-2)
		vel_cor(nx+1) = vel_cor(nx)      
	case default 
		stop 'errorxihua,329'
	end select 

	do ix = 1 , nx + 1
		x_cor(ix) = x(ix) + dt* (vel_cor(ix) + vel(ix))/2.d0 
	enddo 

	do ix = 1 , nx 
		vol_cor(ix) = x_cor(ix+1) - x_cor(ix)
		tau_cor(ix) = vol_cor(ix) / mass(ix)
		den_cor(ix) = 1.d0 / tau_cor(ix)
		ein_cor(ix) = ein(ix) - (pre2(ix)+pre1(ix))*(vol_cor(ix)-vol(ix))/mass(ix)/2.d0
		call den_ein_2_pre_sos(eos(ix),den_cor(ix),ein_cor(ix),pre_cor(ix),sos_cor(ix))
	enddo 

	vol = vol_cor
	den = den_cor
	ein = ein_cor
	pre = pre_cor 
	sos = sos_cor
	x   = x_cor
	vel = vel_cor 

end Subroutine solver_staggered_MJSLSdpdtauDT_details


end module mdu_MJSLS_dpdtau  
