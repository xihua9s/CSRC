!> Last Modified : Sat Feb 25 23:00:50 2023

Module MDU_Initial_kidder
	use common_inicmt
	use common_meshphys
	use fluid_euler_eqn_basic
Contains


	Subroutine Initial_Mesh_Phys_kidder(ini_cmt,mesh_phys)
		implicit none

		type(inicmt)   :: ini_cmt 
		type(meshphys) :: mesh_phys

		integer :: ix,nmax,nu
		real*8  :: r1,r2,pre1,pre2,den1,den2,gamma
		real*8  :: sound1,sound2,tau,ht,ri,rr,den0,u,pre0,den,pre,entropy
		
		stop 'Initial_Mesh_Phys_kidder Not Ready.'
		nmax = ini_cmt%nx

		r1 = 0.9d0
		r2 = 1.d0
		pre1 = 0.1d0
		pre2 = 10.d0
		den2 = 1.d-2

		select case(ini_cmt%cacysp)
		case('Cartesian')
			nu = 1
		case('Cylindrical')
			nu = 2
		case('Spheral')
			nu = 3
		case default
			print*, 'ini_cmt%cacysp=', ini_cmt%cacysp
			stop 'ini_cmt%cacysp is error.'
		end select 
		gamma = 1.d0 + 2.d0/nu

		!> isentropic condition
		den1 = den2 * (pre1/pre2)**(1.d0/gamma)
		print*, 'den1 =',den1

		entropy = pre2 / den2**(gamma)
		print*, 'entropy =', entropy

		!> ---------------------------------------------   
		!> Cheng Juan JCP 14 
		den1    = 6.31d-4
		den2    = 1.d-2   
		Pre1    = 0.1d0
		Pre2    = 10.d0
		entropy = 2.15d4      
		!> ---------------------------------------------   

		sound1 = sqrt(entropy*gamma*den1**(gamma-1.d0))
		sound2 = sqrt(entropy*gamma*den2**(gamma-1.d0))
		tau = sqrt((gamma-1.d0)/2.d0*(r2*r2-r1*r1)/(sound2**2-sound1**2))
		print*, 'tau=',tau

		ht = sqrt(1.d0 - (6.66d-3/tau)**2)

		open(1,file='OutputData/phys_final.tec')   
		write(1,*) 'VARIABLES = "X","Density","Pressure","Velocity","Internal Energy"'
		write(1,*) 'ZONE I=',nmax,',F=POINT'   
		do ix = 0 , nmax
		ri = r1 + (r2 - r1) / nmax * ix
		rr = ht * ri
		den0 = ( (r2*r2-ri*ri)/(r2*r2-r1*r1)*den1**(gamma-1.d0) + &
			(ri*ri-r1*r1)/(r2*r2-r1*r1)*den2**(gamma-1.d0)     )**(1.d0/(gamma-1.d0))
		pre0 = entropy * (den0)**gamma
		den = ht**(-2.d0/(gamma-1.d0)) * den0
		u   = 0.d0
		pre = ht**(-2.d0/(gamma-1.d0)) * pre0
		write(1,909) rr , den , pre , u , pre/den/(gamma-1.d0)
		enddo
		909 format(1x,6D20.8)   


	End Subroutine Initial_Mesh_Phys_kidder


End Module MDU_Initial_kidder
