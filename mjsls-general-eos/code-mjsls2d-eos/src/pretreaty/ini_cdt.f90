

	Subroutine Initial_CMT(ini_cmt)
		use common_inicmt
		implicit none

		type(inicmt) :: ini_cmt

		open(1,file='InputData/case.ini')
		read(1,*) ini_cmt%case_name
		read(1,*) ini_cmt%nx,ini_cmt%ny
		read(1,*) ini_cmt%cfl
		read(1,*) ini_cmt%t_out
		read(1,*) ini_cmt%it_out
		read(1,*) ini_cmt%node_solver
		close(1)

		ini_cmt%it = 0
		ini_cmt%current_time = 0.d0
		ini_cmt%dt = 0.d0

	End Subroutine Initial_CMT

