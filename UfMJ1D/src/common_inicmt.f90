
module common_inicmt
	implicit none

	type inicmt

		Character(99) 	:: case_name
		character(99) 	:: scheme_type
		integer       	:: nx
		real*8  			:: cfl
		real*8  			:: time_stop
		real*8  			:: current_time
		real*8  			:: dt
		integer 			:: it
		integer 			:: it_out
		character(99) 	:: CaCySp
		character(99) 	:: Node_Solver
		character(99) 	:: rp_solver
		logical  		:: lg_output
		integer  		:: n_out

	endtype inicmt

end module common_inicmt
