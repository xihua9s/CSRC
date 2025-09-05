!> Last Modified : 五  3/ 3 22:22:18 2023
module mdu_fluid_euler_eqn_solver_staggered   
use MDU_Initial
use common_inicmt
use common_meshphys
use fluid_euler_eqn_basic
use mdu_Landshoff
use mdu_MJSLS_CC
use mdu_MJSLS_dpdtau
use mdu_MJSLS_Ideal
use mdu_MJSLS_JWL
use mdu_vonN
use mdu_LSDYNA
use mdu_Kurop
use dt
use output
implicit none

contains

Subroutine fluid_euler_eqn_solver_staggered(ini_cmt,mesh_phys)
	implicit none

	type(inicmt)   :: ini_cmt 
	type(meshphys) :: mesh_phys

	select case (ini_cmt%Node_Solver)

	case ('VonN')

		print*, '---> solver_staggered_VonN(ini_cmt,mesh_phys)'
		call solver_staggered_VonN(ini_cmt,mesh_phys)

	case ('Kurop')

		print*, '---> solver_staggered_Kurop(ini_cmt,mesh_phys)'
		call solver_staggered_Kurop(ini_cmt,mesh_phys)

	case ('LSDYNA')

		print*, '---> solver_staggered_lsdyna(ini_cmt,mesh_phys)'
		call solver_staggered_lsdyna(ini_cmt,mesh_phys)

	case ('Landshoff')
		print*, '---> solver_staggered_Landshoff(ini_cmt,mesh_phys)'
		call solver_staggered_Landshoff(ini_cmt,mesh_phys)


	case ('MJSLS-Ideal')

		print*, '---> solver_staggered_MJSLS(ini_cmt,mesh_phys)'
		call solver_staggered_MJSLS_Ideal(ini_cmt,mesh_phys)

	case ('MJSLS-JWL')

		print*, '---> solver_staggered_MJSLS_JWL(ini_cmt,mesh_phys)'
		call solver_staggered_MJSLS_JWL(ini_cmt,mesh_phys)

	case ('MJSLS-CC')

		print*, '---> solver_staggered_MJSLS_CC(ini_cmt,mesh_phys)'
		call solver_staggered_MJSLS_CC(ini_cmt,mesh_phys)

	case ('MJSLSdpdtau')

		print*, '---> solver_staggered_MJSLSdpdtau(ini_cmt,mesh_phys)'
		call solver_staggered_MJSLSdpdtau(ini_cmt,mesh_phys)

	case ('MJSLSdpdtaudt')

		print*, '---> solver_staggered_MJSLSdpdtau(ini_cmt,mesh_phys)'
		call solver_staggered_MJSLSdpdtaudt(ini_cmt,mesh_phys)

	case default

		stop 'Not Ready, xihua25'

	end select

end Subroutine fluid_euler_eqn_solver_staggered   



end module mdu_fluid_euler_eqn_solver_staggered  
