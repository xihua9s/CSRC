!> Last Modified : Sat Feb 25 16:08:50 2023
Module dt
	use common_inicmt
	use common_meshphys
Contains

	Subroutine cal_dt(ini_cmt,mesh_phys)
		implicit none

		type(inicmt)   :: ini_cmt 
		type(meshphys) :: mesh_phys   

		integer  :: ix
		real*8   :: v1,lambda,dt

		dt = 1.d10
		do ix = 1 , ini_cmt%nx
		v1 = mesh_phys%vertex(ix+1) - mesh_phys%vertex(ix)
		lambda = abs(mesh_phys%con_phys(ix,2)/mesh_phys%con_phys(ix,1)) + mesh_phys%sos(ix)
		if( dt > v1/lambda) dt = v1/lambda
		enddo

		ini_cmt%dt = dt * ini_cmt%cfl

		if( ini_cmt%current_time + ini_cmt%dt > ini_cmt%Time_stop ) then
			ini_cmt%dt = ini_cmt%Time_stop - ini_cmt%current_time
		elseif (ini_cmt%current_time + ini_cmt%dt > ini_cmt%Time_stop/ini_cmt%it_out*ini_cmt%n_out) then
			ini_cmt%dt = ini_cmt%Time_stop/ini_cmt%it_out*ini_cmt%n_out - ini_cmt%current_time
			ini_cmt%lg_output = .true.
			ini_cmt%n_out = ini_cmt%n_out + 1
		endif 

		ini_cmt%It = ini_cmt%It + 1
		ini_cmt%current_time = ini_cmt%current_time + ini_cmt%dt
		print*, ini_cmt%It,ini_cmt%dt,ini_cmt%current_time    

	End Subroutine cal_dt  

End Module dt



