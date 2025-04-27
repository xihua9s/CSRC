!> Last Modified : Sat Feb 25 16:08:50 2023
Module MDU_fluid_euler_eqn_solver
	use common_inicmt
	use common_meshphys
	use fluid_euler_eqn_basic
	use mdu_fluid_euler_eqn_solver_staggered
	implicit none

Contains


	Subroutine fluid_euler_eqn_solver(ini_cmt,mesh_phys)
		implicit none

		type(inicmt)   :: ini_cmt 
		type(meshphys) :: mesh_phys
		integer :: ix

		Select Case(ini_cmt%scheme_type)

		Case('centered')

			stop 'centered'

		Case('staggered')

			Call fluid_euler_eqn_solver_staggered(ini_cmt,mesh_phys)

		Case default

			stop 'ini_cmt%scheme_type is Error, pls Check.'

		End Select

	End Subroutine fluid_euler_eqn_solver


End Module MDU_fluid_euler_eqn_solver
