!> Last Modified : Sat Feb 25 16:08:50 2023
Module output
	use common_inicmt
	use common_meshphys
	use fluid_euler_eqn_basic
Contains   

	Subroutine Check_Phys_conservation(ini_cmt,mesh_phys)
		implicit none

		type(inicmt)   :: ini_cmt 
		type(meshphys) :: mesh_phys

		integer :: ix,iy
		real*8  :: momtum,gamma

		momtum = 0.d0

		do ix = 1 , ini_cmt%nx
		momtum = momtum + mesh_phys%mass(ix)*mesh_phys%con_phys(ix,2)/mesh_phys%con_phys(ix,1)
		!call ieos_to_gamma(mesh_phys%eos(ix),gamma)
		do iy = 1 , 3
		if( isnan(mesh_phys%con_phys(ix,iy))) then
			print*, 'mesh_phys%con_phys(ix,iy) is NAN'
			print*, 'ix=',ix,mesh_phys%con_phys(ix,iy), &
				mesh_phys%pre(ix)/mesh_phys%con_phys(ix,1)**gamma
			stop 'NAN, xihua27'
		endif      
		enddo
		enddo

	End Subroutine Check_Phys_conservation   

	Subroutine output_mesh_phys(ini_cmt,mesh_phys)
		implicit none

		type(inicmt)   :: ini_cmt 
		type(meshphys) :: mesh_phys

		integer :: ix
		character(20)   :: fname
		real*8  :: xmid,gamma

		write(fname,'(I20.4)') ini_cmt%n_out
		open(1,file='OutputData/phys_it_'//trim(adjustl(fname))//'.plt')
		write(1,*) 'VARIABLES = "X","Velocity-P","Density","Pressure","Velocity-C","Internal Energy","Entropy"'
		write(1,*) 'ZONE I=',ini_cmt%nx+1,',F=BLOCK,VARLOCATION=([3-7]=CELLCENTERED)'  
		do ix = 1 , ini_cmt%nx + 1
		write(1,*) mesh_phys%vertex(ix)
		enddo 
		do ix = 1 , ini_cmt%nx + 1
		write(1,*) mesh_phys%vel_p(ix)
		enddo 
		do ix = 1 , ini_cmt%nx
		write(1,*) mesh_phys%con_phys(ix,1)
		enddo 
		do ix = 1 , ini_cmt%nx
		write(1,*) mesh_phys%pre(ix)
		enddo 
		do ix = 1 , ini_cmt%nx
		write(1,*) mesh_phys%con_phys(ix,2)/mesh_phys%con_phys(ix,1)
		enddo 
		do ix = 1 , ini_cmt%nx
		write(1,*) mesh_phys%ein(ix)
		enddo 
		do ix = 1 , ini_cmt%nx   
		 !call ieos_to_gamma(mesh_phys%eos(ix),gamma)
		 gamma = 1.4d0
		write(1,*) mesh_phys%pre(ix)/mesh_phys%con_phys(ix,1)**gamma
		enddo
		close(1)

	End Subroutine output_mesh_phys

	Subroutine output_mesh_phys_final(ini_cmt,mesh_phys)
		implicit none

		type(inicmt)   :: ini_cmt 
		type(meshphys) :: mesh_phys

		integer :: ix
		character(8)   :: fname
		real*8  :: xmid,gamma

		print*, '---> output_mesh_phys_final'
		open(1,file='OutputData/phys_final.plt')   
		write(1,'(100A)') 'VARIABLES = "X","Velocity-P","Density","Pressure","Velocity-C","Internal Energy","Entropy"'
		write(1,*) 'ZONE I=',ini_cmt%nx+1,',F=BLOCK,VARLOCATION=([3-7]=CELLCENTERED)'  
		do ix = 1 , ini_cmt%nx + 1
		write(1,*) mesh_phys%vertex(ix)
		enddo 
		do ix = 1 , ini_cmt%nx + 1
		write(1,*) mesh_phys%vel_p(ix)
		enddo 
		do ix = 1 , ini_cmt%nx
		write(1,*) mesh_phys%con_phys(ix,1)
		enddo 
		do ix = 1 , ini_cmt%nx
		write(1,*) mesh_phys%pre(ix)
		enddo 
		do ix = 1 , ini_cmt%nx
		write(1,*) mesh_phys%con_phys(ix,2)/mesh_phys%con_phys(ix,1)
		enddo 
		do ix = 1 , ini_cmt%nx
		write(1,*) mesh_phys%ein(ix)
		enddo 
		do ix = 1 , ini_cmt%nx   
		!call ieos_to_gamma(mesh_phys%eos(ix),gamma)
		gamma = 1.4d0
		write(1,*) mesh_phys%pre(ix)/mesh_phys%con_phys(ix,1)**gamma
		enddo
		close(1)
	End Subroutine output_mesh_phys_final
end Module output
