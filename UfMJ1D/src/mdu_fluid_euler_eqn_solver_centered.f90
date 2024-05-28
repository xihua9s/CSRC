
module mdu_fluid_euler_eqn_solver_centered

	use common_inicmt
	use common_meshphys
	use fluid_euler_eqn_basic
	use MDU_Initial
	use dt
	use output
	use mdu_fluid_euler_eqn_solver_centered_node_solver
	use mdu_fluid_euler_eqn_solver_centered_update
	implicit none

contains


	Subroutine fluid_euler_eqn_solver_centered(ini_cmt,mesh_phys)
		implicit none

		type(inicmt)   :: ini_cmt 
		type(meshphys) :: mesh_phys
		integer :: ix

		mesh_phys%voln_ratio = 0.d0
		Do While (ini_cmt%current_time < ini_cmt%Time_stop)

		call cal_dt(ini_cmt,mesh_phys)

		call fluid_euler_eqn_solver_centered_node_solver(ini_cmt,mesh_phys)

		call update_vertex(ini_cmt,mesh_phys)

		Call update_phys(ini_cmt,mesh_phys)

		if( ini_cmt%lg_output ) then
			Call output_mesh_phys(ini_cmt,mesh_phys)
			ini_cmt%lg_output = .false.
		endif

		Enddo

		Call output_mesh_phys_final(ini_cmt,mesh_phys)

	End Subroutine fluid_euler_eqn_solver_centered

end module mdu_fluid_euler_eqn_solver_centered
