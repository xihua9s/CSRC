!> Last Modified : 五  3/ 3 22:22:18 2023
module mdu_vonN   
use MDU_Initial
use common_inicmt
use common_meshphys
use fluid_euler_eqn_basic
use dt
use output
implicit none

contains



Subroutine solver_staggered_VonN(ini_cmt,mesh_phys)
	implicit none

	type(inicmt)   :: ini_cmt 
	type(meshphys) :: mesh_phys

	ini_cmt%n_out = 1
	ini_cmt%lg_output = .false.
	Do While (ini_cmt%current_time < ini_cmt%Time_stop)

	call cal_dt(ini_cmt,mesh_phys)

	Call solver_staggered_VonN_details( &
		ini_cmt%nx,      &
		ini_cmt%dt,      &
		mesh_phys%con_phys(:,1), &
		mesh_phys%pre, &
		mesh_phys%sos, &
		mesh_phys%ein, &
		mesh_phys%mass, &
		mesh_phys%vol, &
		mesh_phys%vertex, &
		mesh_phys%vel_p &
		)

	if( ini_cmt%lg_output ) then
		Call output_mesh_phys(ini_cmt,mesh_phys)
		ini_cmt%lg_output = .false.
	endif

	Enddo   

	call output_mesh_phys_final(ini_cmt,mesh_phys)

end Subroutine solver_staggered_VonN  


Subroutine solver_staggered_VonN_details(nx,dt,den,pre,sos,ein,mass,vol,x,vel)
	implicit none

	integer,intent(in)   :: nx 
	real*8,intent(in)    :: mass(nx),dt
	real*8,intent(inout) :: den(nx),pre(nx),sos(nx),ein(nx),vol(nx)
	real*8,intent(inout) :: x(nx+1),vel(nx+1)

	integer :: ix
	real*8  :: vis1,vis2,gamma
	real*8  :: q(nx),vel_pre(nx+1),x_pre(nx+1),pre_pre(nx)
	real*8  :: vol_pre(nx),tau_pre(nx),den_pre(nx),ein_pre(nx),sos_pre(nx)
	real*8  :: q_pre(nx),vel_cor(nx+1),x_cor(nx+1),pre_cor(nx)
	real*8  :: vol_cor(nx),tau_cor(nx),den_cor(nx),ein_cor(nx),sos_cor(nx)

	vis1  = 1.d0
	vis2  = 1.d0 
	gamma = 1.4d0 

	do ix = 1 , nx 
	if( vel(ix+1) < vel(ix) ) then
		q(ix) = -vis1*den(ix)*sos(ix)*(vel(ix+1)-vel(ix)) + &
			vis2*den(ix)*(vel(ix+1)-vel(ix))**2
	else
		q(ix) = 0.d0 
	endif 
	enddo 

	vel_pre = vel
	do ix = 2 , nx 
	vel_pre(ix) = vel(ix) - dt*(pre(ix)+q(ix)-pre(ix-1)-q(ix-1)) / (mass(ix)+mass(ix-1))*2.d0
	enddo 

	do ix = 1 , nx + 1
	x_pre(ix) = x(ix) + dt*(vel(ix) + vel_pre(ix))/2.d0
	enddo 

	do ix = 1 , nx 
	vol_pre(ix) = x_pre(ix+1) - x_pre(ix)
	tau_pre(ix) = vol_pre(ix) / mass(ix)
	den_pre(ix) = 1.d0 / tau_pre(ix)
	ein_pre(ix) = ein(ix) - (pre(ix)+q(ix))*(vol_pre(ix)-vol(ix))/mass(ix)
	pre_pre(ix) = (gamma - 1.d0) * den_pre(ix) * ein_pre(ix)
	sos_pre(ix) = sqrt(gamma*pre_pre(ix)/den_pre(ix))
	enddo 

	do ix = 1 , nx 
	if( vel_pre(ix+1) < vel_pre(ix) ) then
		q_pre(ix) = -vis1*den_pre(ix)*sos_pre(ix)*(vel_pre(ix+1)-vel_pre(ix)) + &
			vis2*den_pre(ix)*(vel_pre(ix+1)-vel_pre(ix))**2
	else
		q_pre(ix) = 0.d0 
	endif 
	enddo

	vel_cor = vel_pre
	do ix = 2 , nx 
	vel_cor(ix) = vel(ix) - dt*( pre_pre(ix) + q_pre(ix) - pre_pre(ix-1) - q_pre(ix-1) + &
		pre(ix) + q(ix) - pre(ix-1) - q(ix-1)) / (mass(ix)+mass(ix-1))*2.d0 / 2.d0
	enddo 

	do ix = 1 , nx + 1
	x_cor(ix) = x(ix) + dt* (vel_cor(ix) + vel(ix))/2.d0 
	enddo 

	do ix = 1 , nx 
	vol_cor(ix) = x_cor(ix+1) - x_cor(ix)
	tau_cor(ix) = vol_cor(ix) / mass(ix)
	den_cor(ix) = 1.d0 / tau_cor(ix)
	ein_cor(ix) = ein(ix) - (pre_pre(ix)+q_pre(ix)+pre(ix)+q(ix))*(vol_cor(ix)-vol(ix))/mass(ix)/2.d0
	pre_cor(ix) = (gamma - 1.d0)*den_cor(ix)*ein_cor(ix)
	sos_cor(ix) = sqrt(gamma*pre_cor(ix)/den_cor(ix))
	enddo 

	vol = vol_cor
	den = den_cor
	ein = ein_cor
	pre = pre_cor 
	sos = sos_cor
	x   = x_pre
	vel = vel_cor 


end Subroutine solver_staggered_VonN_details


end module mdu_vonN  
