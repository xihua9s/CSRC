
Module fluid_euler_eqn_basic   
	use mdu_ideal 
	use mdu_cc 
	use mdu_jwl
	implicit none

contains   


subroutine PVarToSoS( nUvar, PVar,sos, iEOS)
	implicit none  
	real*8, intent(out):: sos
	real*8, dimension(nUvar),intent(in)::PVar
	real*8  :: gamma6(6),d,p,ftau,ein,Gamma
	integer :: iEOS,nUvar,n
	if(iEOS < 0 ) stop 'iEOS is Error'   

	Select case (ieos)
	case(140)
		gamma = 1.4d0
		call PVarToOther( nUvar, PVar,sos,ein, Gamma , iEOS)
	case(166)
		gamma = 5.d0/3.d0
		call PVarToOther( nUvar, PVar,sos,ein, Gamma , iEOS)
	case(715)
		gamma = 7.d0 / 5.d0
		call PVarToOther( nUvar, PVar,sos,ein, Gamma , iEOS)
	case(513)
		gamma = 5.d0 / 3.d0   
		call PVarToOther( nUvar, PVar,sos,ein, Gamma , iEOS)
	case(17)					
		call PVarToOther_JWL( nUvar, PVar,sos,ein, iEOS)
	case(9502)					
		call PVarToOther_JWL( nUvar, PVar,sos,ein, iEOS)
	case(95021)					
		call PVarToOther_JWL( nUvar, PVar,sos,ein, iEOS)
	case(95022)					
		call PVarToOther_JWL( nUvar, PVar,sos,ein, iEOS)
	case(33)					
		call PVarToOther_CC( nUvar, PVar,sos,ein, iEOS)
	case default
		stop '32fa,xihua'
	End select 

	if ( sos <= 0.d0 ) then 
		print*, sos ,ein
		stop 'sos is less than 0bb'
	else 
		sos = sqrt(sos)
	endif  
endsubroutine PVarToSoS


subroutine den_ein_2_pre_sos(iEOS,d,ein,p,sos)
	implicit none
	integer :: iEOS 
	real*8,intent(in)  :: d,ein
	real*8,intent(out) :: p,sos
	real*8 :: gamma6(6),gamma

	Select case (ieos)
	case(140)
		gamma = 1.4d0
		p = (gamma - 1.d0) * d * ein
		sos = sqrt(gamma*p/d)
	case(166)
		gamma = 5.d0/3.d0
		p = (gamma - 1.d0) * d * ein
		sos = sqrt(gamma*p/d)
	case(715)
		gamma = 7.d0 / 5.d0
		p = (gamma - 1.d0) * d * ein
		sos = sqrt(gamma*p/d)
	case(513)
		gamma = 5.d0 / 3.d0   
		p = (gamma - 1.d0) * d * ein
		sos = sqrt(gamma*p/d)
	case(17)					
		call den_ein_2_pre_sos_JWL(iEOS,d,ein,p,sos)		
	case(9502)					
		call den_ein_2_pre_sos_JWL(iEOS,d,ein,p,sos)		
	case(95021)					
		call den_ein_2_pre_sos_JWL(iEOS,d,ein,p,sos)		
	case(95022)					
		call den_ein_2_pre_sos_JWL(iEOS,d,ein,p,sos)		
	case(33)			
		call den_ein_2_pre_sos_CC(iEOS,d,ein,p,sos)		
	case default
		stop '32fa,xihua'
	End select 

end subroutine den_ein_2_pre_sos

subroutine ieos_to_dpdtau(ieos,tau,pre,dpdtau,dp2dtau2)
	implicit none 

	integer,intent(in)	:: ieos 
	real*8,intent(in)		:: tau,pre 
	real*8,intent(out)	:: dpdtau,dp2dtau2

	real*8	:: gamma,gamma6(6)
	select case(iEos)
	case(140)
		gamma = 1.4d0
		dpdtau = -gamma*pre/tau
		dp2dtau2 = gamma*(gamma+1.d0)/2.d0*pre/tau/tau
	case(166)
		gamma = 5.d0/3.d0
		dpdtau = -gamma*pre/tau
		dp2dtau2 = gamma*(gamma+1.d0)/2.d0*pre/tau/tau
	case(715)
		gamma = 7.d0 / 5.d0
		dpdtau = -gamma*pre/tau
		dp2dtau2 = gamma*(gamma+1.d0)/2.d0*pre/tau/tau
	case(513)
		gamma = 5.d0 / 3.d0     
		dpdtau = -gamma*pre/tau
		dp2dtau2 = gamma*(gamma+1.d0)/2.d0*pre/tau/tau
	case(17) !LX17
		call ieos_to_gamma_JWL(ieos,gamma6)
		dpdtau = (gamma6(3)*(1.d0+gamma6(1)-gamma6(5)*gamma6(2)*tau) * &
					 exp(-gamma6(5)*gamma6(2)*tau) + &
		          gamma6(4)*(1.d0+gamma6(1)-gamma6(6)*gamma6(2)*tau) * &
					 exp(-gamma6(6)*gamma6(2)*tau) - &
					 (gamma6(1)+1.d0)*pre)/tau
		dp2dtau2 =-(2.d0+gamma6(1))*(gamma6(3)*(1+gamma6(1)- &
						gamma6(5)*gamma6(2)*tau) &
						*exp(-gamma6(5)*gamma6(2)*tau) + &
						gamma6(4)*(1+gamma6(1)- &
						gamma6(6)*gamma6(2)*tau) &
						*exp(-gamma6(6)*gamma6(2)*tau))/tau/tau &
						-(gamma6(3)*gamma6(5)*gamma6(2)*(2.d0+gamma6(1) - &
						gamma6(5)*gamma6(2)*tau) &
						*exp(-gamma6(5)*gamma6(2)*tau) &
						+gamma6(4)*gamma6(6)*gamma6(2)*(2.d0+gamma6(1) - &
						gamma6(6)*gamma6(2)*tau) &
						*exp(-gamma6(6)*gamma6(2)*tau))/tau &
						+(2.d0+gamma6(1))*(1.d0+gamma6(1))*pre/tau/tau		
	case(33) !CC-33
		call ieos_to_gamma_CC(ieos,gamma6)
		dpdtau = ( gamma6(3)*(1.d0-gamma6(5)+gamma6(1)) * &
		         (gamma6(2)*tau)**(-gamma6(5)) &
					-gamma6(4)*(1.d0-gamma6(6)+gamma6(1)) * &
					(gamma6(2)*tau)**(-gamma6(6)) &
					-(1.d0+gamma6(1))*pre ) / tau
		dp2dtau2 = ((2.d0+gamma6(1))*(1.d0+gamma6(1))*pre - &
						gamma6(3)*(2.d0+gamma6(5)+gamma6(1))* &
						(1.d0-gamma6(5)+gamma6(1))* &
						(gamma6(2)*tau)**(-gamma6(5)))/tau/tau &
					 + gamma6(4)*(2.d0+gamma6(6)+gamma6(1)) * &
					 (1.d0-gamma6(6)+gamma6(1))* &
					 (gamma6(2)*tau)**(-gamma6(6))/tau/tau

	case(9502)  !PBX9502
		call ieos_to_gamma_JWL(ieos,gamma6)
		dpdtau = (gamma6(3)*(1.d0+gamma6(1)-gamma6(5)*gamma6(2)*tau) * &
					 exp(-gamma6(5)*gamma6(2)*tau) + &
		          gamma6(4)*(1.d0+gamma6(1)-gamma6(6)*gamma6(2)*tau) * &
					 exp(-gamma6(6)*gamma6(2)*tau) - &
					 (gamma6(1)+1.d0)*pre)/tau
		dp2dtau2 =-(2.d0+gamma6(1))*(gamma6(3)*(1+gamma6(1)- &
						gamma6(5)*gamma6(2)*tau) &
						*exp(-gamma6(5)*gamma6(2)*tau) + &
						gamma6(4)*(1+gamma6(1)- &
						gamma6(6)*gamma6(2)*tau) &
						*exp(-gamma6(6)*gamma6(2)*tau))/tau/tau &
						-(gamma6(3)*gamma6(5)*gamma6(2)*(2.d0+gamma6(1) - &
						gamma6(5)*gamma6(2)*tau) &
						*exp(-gamma6(5)*gamma6(2)*tau) &
						+gamma6(4)*gamma6(6)*gamma6(2)*(2.d0+gamma6(1) - &
						gamma6(6)*gamma6(2)*tau) &
						*exp(-gamma6(6)*gamma6(2)*tau))/tau &
						+(2.d0+gamma6(1))*(1.d0+gamma6(1))*pre/tau/tau		
	case(95021) !PBX9502_Sod
		call ieos_to_gamma_JWL(ieos,gamma6)
		dpdtau = (gamma6(3)*(1.d0+gamma6(1)-gamma6(5)*gamma6(2)*tau) * &
					 exp(-gamma6(5)*gamma6(2)*tau) + &
		          gamma6(4)*(1.d0+gamma6(1)-gamma6(6)*gamma6(2)*tau) * &
					 exp(-gamma6(6)*gamma6(2)*tau) - &
					 (gamma6(1)+1.d0)*pre)/tau
		dp2dtau2 =-(2.d0+gamma6(1))*(gamma6(3)*(1+gamma6(1)- &
						gamma6(5)*gamma6(2)*tau) &
						*exp(-gamma6(5)*gamma6(2)*tau) + &
						gamma6(4)*(1+gamma6(1)- &
						gamma6(6)*gamma6(2)*tau) &
						*exp(-gamma6(6)*gamma6(2)*tau))/tau/tau &
						-(gamma6(3)*gamma6(5)*gamma6(2)*(2.d0+gamma6(1) - &
						gamma6(5)*gamma6(2)*tau) &
						*exp(-gamma6(5)*gamma6(2)*tau) &
						+gamma6(4)*gamma6(6)*gamma6(2)*(2.d0+gamma6(1) - &
						gamma6(6)*gamma6(2)*tau) &
						*exp(-gamma6(6)*gamma6(2)*tau))/tau &
						+(2.d0+gamma6(1))*(1.d0+gamma6(1))*pre/tau/tau		
	case(95022) !PBX9502_WC
		call ieos_to_gamma_JWL(ieos,gamma6)
		dpdtau = (gamma6(3)*(1.d0+gamma6(1)-gamma6(5)*gamma6(2)*tau) * &
					 exp(-gamma6(5)*gamma6(2)*tau) + &
		          gamma6(4)*(1.d0+gamma6(1)-gamma6(6)*gamma6(2)*tau) * &
					 exp(-gamma6(6)*gamma6(2)*tau) - &
					 (gamma6(1)+1.d0)*pre)/tau
		dp2dtau2 =-(2.d0+gamma6(1))*(gamma6(3)*(1+gamma6(1)- &
						gamma6(5)*gamma6(2)*tau) &
						*exp(-gamma6(5)*gamma6(2)*tau) + &
						gamma6(4)*(1+gamma6(1)- &
						gamma6(6)*gamma6(2)*tau) &
						*exp(-gamma6(6)*gamma6(2)*tau))/tau/tau &
						-(gamma6(3)*gamma6(5)*gamma6(2)*(2.d0+gamma6(1) - &
						gamma6(5)*gamma6(2)*tau) &
						*exp(-gamma6(5)*gamma6(2)*tau) &
						+gamma6(4)*gamma6(6)*gamma6(2)*(2.d0+gamma6(1) - &
						gamma6(6)*gamma6(2)*tau) &
						*exp(-gamma6(6)*gamma6(2)*tau))/tau &
						+(2.d0+gamma6(1))*(1.d0+gamma6(1))*pre/tau/tau		
	case default
		stop 'ieos is error, line11'
	end select

end subroutine ieos_to_dpdtau


End Module fluid_euler_eqn_basic
