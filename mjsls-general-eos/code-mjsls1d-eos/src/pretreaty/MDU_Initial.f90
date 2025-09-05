!> Last Modified : Sat Feb 25 22:30:15 2023

Module MDU_Initial
	use common_inicmt
	use common_meshphys
	use fluid_euler_eqn_basic
	use MDU_Initial_sod
	use MDU_Initial_lax
	use MDU_Initial_shuosher
	use MDU_Initial_twoblast
	use MDU_Initial_leblanc
	use MDU_Initial_sedov
	use MDU_Initial_noh
	use MDU_Initial_s123
	use MDU_Initial_s321
	use MDU_Initial_lazarus
	use MDU_Initial_kidder
	use MDU_Initial_JWL
	use MDU_Initial_stiff
	use MDU_Initial_TNT
Contains

	Subroutine  Initial(ini_cmt,mesh_phys)
		implicit none

		type(inicmt)   :: ini_cmt 
		type(meshphys) :: mesh_phys

		Call Initial_CMT(ini_cmt)

		Call rege_mesh_phys(ini_cmt,mesh_phys)

		Call Initial_Mesh_Phys(ini_cmt,mesh_phys)

	End Subroutine  Initial

	Subroutine Initial_CMT(ini_cmt)
		implicit none

		type(inicmt) :: ini_cmt

		open(1,file='InputData/case.ini')
		read(1,*) ini_cmt%case_name
		read(1,*) ini_cmt%scheme_type
		read(1,*) ini_cmt%CaCySp
		read(1,*) ini_cmt%nx
		read(1,*) ini_cmt%cfl
		read(1,*) ini_cmt%it_out
		read(1,*) ini_cmt%Node_Solver
		close(1)

		ini_cmt%it = 0
		ini_cmt%current_time = 0.d0
		ini_cmt%dt = 0.d0

	End Subroutine Initial_CMT

	Subroutine rege_mesh_phys(ini_cmt,mesh_phys)
		implicit none

		type(inicmt)   :: ini_cmt 
		type(meshphys) :: mesh_phys

		Select Case(ini_cmt%scheme_type)

		Case('centered')

			allocate(mesh_phys%vertex(ini_cmt%nx+1))
			allocate(mesh_phys%vel_p(ini_cmt%nx+1))
			allocate(mesh_phys%pstar(ini_cmt%nx+1))

			allocate(mesh_phys%con_phys(ini_cmt%nx,3))
			allocate(mesh_phys%pre(ini_cmt%nx))
			allocate(mesh_phys%sos(ini_cmt%nx))
			allocate(mesh_phys%ein(ini_cmt%nx))
			allocate(mesh_phys%eos(ini_cmt%nx))
			allocate(mesh_phys%mass(ini_cmt%nx))
			allocate(mesh_phys%vol_size_t0(ini_cmt%nx))
			allocate(mesh_phys%vol_size_old(ini_cmt%nx))
			allocate(mesh_phys%vol_size_new(ini_cmt%nx))
			allocate(mesh_phys%entropy(ini_cmt%nx))
			allocate(mesh_phys%volnmin1(ini_cmt%nx))
			allocate(mesh_phys%vol0(ini_cmt%nx))
			allocate(mesh_phys%den0(ini_cmt%nx))
			allocate(mesh_phys%voln_ratio(ini_cmt%nx))
			allocate(mesh_phys%tau(ini_cmt%nx))
			allocate(mesh_phys%tau_old(ini_cmt%nx))  
			allocate(mesh_phys%vol(ini_cmt%nx))          

			mesh_phys%entropy = 0.d0

		Case('staggered')

			allocate(mesh_phys%vertex(ini_cmt%nx+1))
			allocate(mesh_phys%vel_p(ini_cmt%nx+1))

			allocate(mesh_phys%pstar(ini_cmt%nx))   
			allocate(mesh_phys%con_phys(ini_cmt%nx,3))
			allocate(mesh_phys%pre(ini_cmt%nx))
			allocate(mesh_phys%sos(ini_cmt%nx))
			allocate(mesh_phys%ein(ini_cmt%nx))
			allocate(mesh_phys%eos(ini_cmt%nx))
			allocate(mesh_phys%mass(ini_cmt%nx))
			allocate(mesh_phys%vol_size_old(ini_cmt%nx))
			allocate(mesh_phys%vol_size_t0(ini_cmt%nx))
			allocate(mesh_phys%vol(ini_cmt%nx))
			allocate(mesh_phys%vol_size_new(ini_cmt%nx))
			allocate(mesh_phys%tau(ini_cmt%nx))
			allocate(mesh_phys%tau_old(ini_cmt%nx))
			allocate(mesh_phys%q(ini_cmt%nx))

		Case default

			stop 'ini_cmt%scheme_type is Error, pls Check.'

		End Select      

	end Subroutine rege_mesh_phys

	Subroutine Initial_Mesh_Phys(ini_cmt,mesh_phys)
		implicit none

		type(inicmt)   :: ini_cmt 
		type(meshphys) :: mesh_phys
		integer  :: ix,nu

		Select Case(ini_cmt%case_name)

		Case('Sod')

			Call Initial_Mesh_Phys_sod(ini_cmt,mesh_phys)

		Case('Lax')

			Call Initial_Mesh_Phys_lax(ini_cmt,mesh_phys)

		Case('ShuOsher')

			Call Initial_Mesh_Phys_shuosher(ini_cmt,mesh_phys)

		Case('TwoBlast')

			Call Initial_Mesh_Phys_twoblast(ini_cmt,mesh_phys)

		Case('Leblanc')

			Call Initial_Mesh_Phys_Leblanc(ini_cmt,mesh_phys)      

		Case('S123')

			Call Initial_Mesh_Phys_s123(ini_cmt,mesh_phys)

		Case('S321')

			Call Initial_Mesh_Phys_s321(ini_cmt,mesh_phys)

		Case('Noh')

			Call Initial_Mesh_Phys_noh(ini_cmt,mesh_phys)

		Case('Sedov')

			Call Initial_Mesh_Phys_sedov(ini_cmt,mesh_phys)

		Case('Lazarus')

			Call Initial_Mesh_Phys_Lazarus(ini_cmt,mesh_phys)   

		Case('Kidder')

			Call Initial_Mesh_Phys_Kidder(ini_cmt,mesh_phys)         

		Case('Exp2vac')

			Call Initial_Mesh_Phys_Exp2vac(ini_cmt,mesh_phys) 

		Case('Water-air')

			Call Initial_Mesh_Phys_water_air(ini_cmt,mesh_phys) 

		Case('JWL-LX17')

			Call Initial_Mesh_Phys_JWL_LX17(ini_cmt,mesh_phys) 

		Case('JWL-LX17-Blast')

			Call Initial_Mesh_Phys_JWL_LX17_Blast(ini_cmt,mesh_phys) 

		Case('JWL-PBX9502')

			Call Initial_Mesh_Phys_JWL_PBX9502(ini_cmt,mesh_phys) 

		Case('JWL-PBX9502-Sod')

			Call Initial_Mesh_Phys_JWL_PBX9502_Sod(ini_cmt,mesh_phys) 

		Case('TNT-water')

			Call Initial_Mesh_Phys_TNT_water(ini_cmt,mesh_phys) 

		Case('CC-33')

			Call Initial_Mesh_Phys_CC_33(ini_cmt,mesh_phys) 

		Case('vandeW-Sod')

			Call Initial_Mesh_Phys_vandeW_Sod(ini_cmt,mesh_phys) 

		Case('Wilkins')
print*, 1
			Call Initial_Mesh_Phys_Wilkins(ini_cmt,mesh_phys) 
print*, 2
		Case default

			print*, 'ini_cmt%case_name=',ini_cmt%case_name
			stop 'ini_cmt%case_name is error, pls call Xihua'

		End Select

		select case(ini_cmt%cacysp)
		case('Cartesian')
			nu = 0
		case('Cylindrical')
			nu = 1
		case('Spheral')
			nu = 2
		case default
			print*, 'ini_cmt%cacysp=', ini_cmt%cacysp
			stop 'ini_cmt%cacysp is error.'
		end select    

		do ix = 1 , ini_cmt%nx      
			mesh_phys%vol(ix) = mesh_phys%vertex(ix+1)**(nu+1)/(nu+1) - &
				mesh_phys%vertex(ix)**(nu+1)/(nu+1)    
			mesh_phys%mass(ix) = mesh_phys%con_phys(ix,1) * mesh_phys%vol(ix)   
			mesh_phys%vol_size_t0(ix)  = mesh_phys%vol_size_old(ix)    
			mesh_phys%tau(ix) = 1.d0 / mesh_phys%con_phys(ix,1)    
			mesh_phys%tau_old(ix) =  mesh_phys%tau(ix)  
		enddo   

	End Subroutine Initial_Mesh_Phys


End Module MDU_Initial
