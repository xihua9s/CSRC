
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
	case(166,513)
		gamma = 5.d0/3.d0
		call PVarToOther( nUvar, PVar,sos,ein, Gamma , iEOS)
	case(715)
		gamma = 7.d0 / 5.d0
		call PVarToOther( nUvar, PVar,sos,ein, Gamma , iEOS)
	case(17,9502,95021,95022,3712)					
		call PVarToOther_JWL( nUvar, PVar,sos,ein, iEOS)
	case(33)					
		call PVarToOther_CC( nUvar, PVar,sos,ein, iEOS)
	case(1920,7153309) !> Stiffened gas EOS
		call PvarToOther_Stiff( nUvar, PVar,sos,ein, iEOS)
	case(22423)
		call PVarToOther_vandeW( nUvar, PVar,sos,ein, iEOS)
	case(14719)
		call PVarToOther_MieGruSolid( nUvar, PVar,sos,ein, iEOS)
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
	case(166,513)
		gamma = 5.d0/3.d0
		p = (gamma - 1.d0) * d * ein
		sos = sqrt(gamma*p/d)
	case(715)
		gamma = 7.d0 / 5.d0
		p = (gamma - 1.d0) * d * ein
		sos = sqrt(gamma*p/d)
	case(17,9502,95021,95022,3712)					
		call den_ein_2_pre_sos_JWL(iEOS,d,ein,p,sos)		
	case(33)			
		call den_ein_2_pre_sos_CC(iEOS,d,ein,p,sos)	
	case(1920,7153309)			
		call den_ein_2_pre_sos_stiff(iEOS,d,ein,p,sos)	
	case(22423)
		call den_ein_2_pre_sos_vandeW(iEOS,d,ein,p,sos)	
	case(14719)
		call den_ein_2_pre_sos_MieGruSolid(iEOS,d,ein,p,sos)	
	case default
		stop '32fa,xihua'
	End select 

end subroutine den_ein_2_pre_sos

subroutine ieos_to_dpdtau(ieos,tau,pre,dpdtau,dp2dtau2)
	implicit none 

	integer,intent(in)	:: ieos 
	real*8,intent(in)		:: tau,pre 
	real*8,intent(out)	:: dpdtau,dp2dtau2

	real*8	:: gamma,gamma6(6),gamma3(3),gamma4(4)
	real*8  :: alpha,beta,R,Cv,rho_0,a_0,gma_0,Sm,tau_0
	select case(iEos)
	case(140)
		gamma = 1.4d0
		dpdtau = -gamma*pre/tau
		dp2dtau2 = gamma*(gamma+1.d0)*pre/tau/tau
	case(166,513)
		gamma = 5.d0/3.d0
		dpdtau = -gamma*pre/tau
		dp2dtau2 = gamma*(gamma+1.d0)*pre/tau/tau
	case(715)
		gamma = 7.d0 / 5.d0
		dpdtau = -gamma*pre/tau
		dp2dtau2 = gamma*(gamma+1.d0)*pre/tau/tau
	case(17,9502,95021,95022,3712) !JWL
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
	case (1920,7153309) !> stiff
		call ieos_to_gamma_stiff(ieos,gamma3)
		dpdtau = -gamma3(1)*(pre+gamma3(3))/tau 
		dp2dtau2 = gamma3(1)*(gamma3(1)+1.d0)*(pre+gamma3(3))/tau/tau
	case(22423) !> van der Walls
		call ieos_to_gamma_vandeW(ieos,gamma4)
		alpha = gamma4(1)
		beta  = gamma4(2)
		R     = gamma4(3)
		Cv    = gamma4(4)
		dpdtau = (2.d0*Cv*alpha*(1.d0-beta/tau) - &
				 (Cv+R)*(tau*tau*pre+alpha)) / &
				 (Cv*tau*tau*tau-Cv*beta*tau*tau)
		dp2dtau2 = 2.d0*( Cv*(tau-beta)*pre + &
		                  Cv*alpha*(tau-beta)/tau/tau - &
					R*alpha/tau - R*tau*pre ) / &
					(Cv*tau**3-Cv*beta*tau*tau) &
				-0.5d0*(4.d0*(Cv*tau*(tau-beta)*pre+ &
				Cv*alpha*(1.d0-beta/tau))-2.d0*R*alpha &
				-2.d0*Cv*alpha-R*tau*tau*pre ) / &
				(Cv*tau**3-Cv*beta*tau*tau)**2 * &
				(6.d0*Cv*tau*tau-4.d0*Cv*beta*tau+R*tau*tau) &
				-0.5d0*(Cv*tau*tau*(tau-beta)*pre)/ &
				(Cv*tau**3-Cv*beta*tau*tau)**2* &
				(12.d0*Cv*tau-4.d0*Cv*beta+4.d0*R*tau) + &
				0.5d0*(Cv*tau*tau*(tau-beta)*pre)/ &
				(Cv*tau**3-Cv*beta*tau*tau)**3* &
				(6.d0*Cv*tau*tau-4.d0*Cv*beta*tau+R*tau*tau)**2
	case(14719) !> Mie-Gruneisen EOS for solids
		Call ieos_to_gamma_MieGruSolid(ieos,gamma4)
		rho_0 = gamma4(1)
		a_0   = gamma4(2)
		gma_0 = gamma4(3)
		Sm    = gamma4(4)
		tau_0 = 1.d0/rho_0 
		dpdtau = -gma_0*pre/tau_0+a_0*a_0/tau_0*(gma_0*(tau_0-tau)-tau_0)/ &
				 (tau_0-sm*(tau_0-tau))**2-2.d0*a_0*a_0*sm/tau_0*  &
				(tau_0-tau)*(tau_0-0.5d0*gma_0*(tau_0-tau))/(tau_0-sm*(tau_0-tau))**3
		dp2dtau2 = gma_0**2*pre/tau_0/tau_0 - (a_0*a_0*gma_0/tau_0+ &
					a_0*a_0*gma_0/tau_0/tau_0*(gma_0*(tau_0-tau)-tau_0)) / &
					(tau_0-sm*(tau_0-tau))**2 + (2.d0*a_0*a_0*sm*gma_0/tau_0/tau_0* &
					(tau_0-tau)*(tau_0-0.5d0*gma_0*(tau_0-tau))-4.d0*a_0*a_0*sm/tau_0* &
					(gma_0*(tau_0-tau)-tau_0))/(tau_0-sm*(tau_0-tau))**3 + &
					6.d0*a_0*a_0*sm*sm/tau_0*(tau_0-tau)*(tau_0-0.5d0*gma_0*(tau_0-tau))/ &
					(tau_0 - sm*(tau_0-tau))**4
	case default
		stop 'ieos is error, line11'
	end select

end subroutine ieos_to_dpdtau


End Module fluid_euler_eqn_basic
