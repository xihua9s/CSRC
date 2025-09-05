
Module mdu_cc   
	implicit none

contains   

subroutine PVarToOther_CC( nUvar, PVar,sos,ein, iEOS)
	implicit none  
	real*8, intent(out):: sos,ein
	real*8, dimension(nUvar),intent(in)::PVar
	real*8  :: gamma6(6),d,p,ftau
	integer :: iEOS,nUvar,n
	if(iEOS < 0 ) stop 'iEOS is Error'   

	Call ieos_to_gamma_CC(ieos,gamma6)
	n = size(PVar)
	d = PVar(1)
	p = PVar(n)
	ftau = gamma6(3)*(1.d0-gamma6(5)+gamma6(1))/(1.d0-gamma6(5)) * &
			(gamma6(2)/d)**(-gamma6(5)) - &
				gamma6(4)*(1.d0-gamma6(6)+gamma6(1))/(1.d0-gamma6(6)) * &
			(gamma6(2)/d)**(-gamma6(6)) 
	ein = ( p - ftau )/gamma6(1)/d
	sos = (gamma6(1)+1.d0)*p/d - &
			gamma6(3)/d*(1.d0-gamma6(5)+gamma6(1))*(gamma6(2)/d)**(-gamma6(5)) &
			+gamma6(4)/d*(1.d0-gamma6(6)+gamma6(1))*(gamma6(2)/d)**(-gamma6(6))	

	if ( sos <= 0.d0 ) then 
		print*, sos ,ein
		stop 'sos is less than 0bb'
	else 
		sos = sqrt(sos)
	endif  

endsubroutine PVarToOther_CC


subroutine den_ein_2_pre_sos_CC(iEOS,d,ein,p,sos)
	implicit none
	integer :: iEOS 
	real*8,intent(in)  :: d,ein
	real*8,intent(out) :: p,sos
	real*8 :: gamma6(6),ftau 

	Call ieos_to_gamma_CC(ieos,gamma6)

	ftau = gamma6(3)*(1.d0-gamma6(5)+gamma6(1))/(1.d0-gamma6(5)) * &
			(gamma6(2)/d)**(-gamma6(5)) - &
				gamma6(4)*(1.d0-gamma6(6)+gamma6(1))/(1.d0-gamma6(6)) * &
			(gamma6(2)/d)**(-gamma6(6))
	p = gamma6(1)*ein*d + ftau 
	sos = (gamma6(1)+1.d0)*p/d - &
			gamma6(3)/d*(1.d0-gamma6(5)+gamma6(1))*(gamma6(2)/d)**(-gamma6(5)) &
			+gamma6(4)/d*(1.d0-gamma6(6)+gamma6(1))*(gamma6(2)/d)**(-gamma6(6))	
			
	if ( sos <= 0.d0 ) then 
		print*, sos ,ein
		stop 'sos is less than 0bb'
	else 
		sos = sqrt(sos)
	endif  

end subroutine den_ein_2_pre_sos_CC

subroutine CVarToPVar_CC(nuvar,CVar, PVar , iEOS)
	implicit none
	real*8, dimension(nuvar),intent(in)::CVar
	real*8, dimension(nuvar),intent(inout)::PVar
	real*8  :: gamma6(6),d,p,e,ein,ftau 
	integer :: iEOS,n,nuvar
	if(iEOS < 0 ) stop 'iEOS is Error'

	Call ieos_to_gamma_CC(ieos,gamma6)

	n = size(CVar)
	d = CVar(1)
	E = CVar(n)   
	ein = (E-0.5*dot_product(CVar(2:n-1),CVar(2:n-1))/d) / d

	ftau = gamma6(3)*(1.d0-gamma6(5)+gamma6(1))/(1.d0-gamma6(5)) * &
			(gamma6(2)/d)**(-gamma6(5)) - &
				gamma6(4)*(1.d0-gamma6(6)+gamma6(1))/(1.d0-gamma6(6)) * &
			(gamma6(2)/d)**(-gamma6(6))
	p = gamma6(1)*ein*d + ftau 

	PVar(1) = d
	PVar(2:n-1) = CVar(2:n-1)/d
	PVar(n) = p

endsubroutine CVarToPVar_CC


subroutine PVarToCVar_CC( nUvar, PVar,CVar , iEOS)
	implicit none  
	real*8, dimension(nUvar),intent(inout)::CVar
	real*8, dimension(nUvar),intent(in)::PVar
	real*8  :: gamma6(6),d,p,e,ein,ftau
	integer :: iEOS,n,nUvar
	if(iEOS < 0 ) stop 'iEOS is Error'   

	Call ieos_to_gamma_CC(ieos,gamma6)
	n = size(PVar)
	d = PVar(1)
	p = PVar(n)

	ftau = gamma6(3)*(1.d0-gamma6(5)+gamma6(1))/(1.d0-gamma6(5)) * &
			(gamma6(2)/d)**(-gamma6(5)) - &
				gamma6(4)*(1.d0-gamma6(6)+gamma6(1))/(1.d0-gamma6(6)) * &
			(gamma6(2)/d)**(-gamma6(6)) 
	ein = ( p - ftau )/gamma6(1)/d

	E = 0.5*d*dot_product(PVar(2:n-1),PVar(2:n-1)) + ein*d
	CVar(1) = d
	CVar(2:n-1) = PVar(2:n-1)*d
	CVar(n) = E
endsubroutine PVarToCVar_CC

subroutine ieos_to_gamma_CC(ieos,gamma6)
	implicit none 
	integer,intent(in)   :: ieos 
	real*8,intent(out)   :: gamma6(6)

	gamma6 = 0.d0 
	Select case(ieos)
	case(33)
		gamma6(1) = 1.19d0
		gamma6(2) = 1134.d0 
		gamma6(3) = 8.192d9
		gamma6(4) = 1.508d9
		gamma6(5) = 4.53d0 
		gamma6(6) = 1.42d0 
	case default
		stop 'CC EOS is error.'
	End Select
end subroutine ieos_to_gamma_CC




subroutine PVarToOther_vandeW( nUvar, PVar,sos,ein, iEOS)
	implicit none  
	real*8, intent(out):: sos,ein
	real*8, dimension(nUvar),intent(in)::PVar
	real*8  :: gamma4(4),d,p,ftau,tau,gg
	integer :: iEOS,nUvar,n
	if(iEOS < 0 ) stop 'iEOS is Error'   

	Call ieos_to_gamma_vandeW(ieos,gamma4)
	n = size(PVar)
	d = PVar(1)
	tau = 1/d
	p = PVar(n)

	ftau = gamma4(3)/gamma4(4)*(gamma4(1)/tau/(tau-gamma4(2)))-gamma4(1)/tau/tau 
	gg   = gamma4(3)/gamma4(4)*(tau/(tau-gamma4(2)))
	ein = ( p - ftau )/gg*tau 

	sos = gamma4(3)/gamma4(4)*((p*tau*tau+ein*tau+2.d0*gamma4(1))/(tau-gamma4(2)) &
			+ (gamma4(2)*ein*tau+gamma4(1)*gamma4(2))/(tau-gamma4(2))**2 )        &
			- 2.d0*gamma4(1)/tau

	if ( sos <= 0.d0 ) then 
		print*, sos ,ein
		stop 'sos is less than 0bb'
	else 
		sos = sqrt(sos)
	endif  

endsubroutine PVarToOther_vandeW


subroutine den_ein_2_pre_sos_vandeW(iEOS,d,ein,p,sos)
	implicit none
	integer :: iEOS 
	real*8,intent(in)  :: d,ein
	real*8,intent(out) :: p,sos
	real*8 :: gamma4(4),ftau,gg,tau

	Call ieos_to_gamma_vandeW(ieos,gamma4)

	tau = 1/d 
	ftau = gamma4(3)/gamma4(4)*(gamma4(1)/tau/(tau-gamma4(2)))-gamma4(1)/tau/tau 
	gg   = gamma4(3)/gamma4(4)*(tau/(tau-gamma4(2)))

	p = gg * ein * d + ftau 
	sos = gamma4(3)/gamma4(4)*((p*tau*tau+ein*tau+2.d0*gamma4(1))/(tau-gamma4(2)) &
			+ (gamma4(2)*ein*tau+gamma4(1)*gamma4(2))/(tau-gamma4(2))**2 )        &
			- 2.d0*gamma4(1)/tau

	if ( sos <= 0.d0 ) then 
		print*, sos ,ein
		stop 'sos is less than 0bb'
	else 
		sos = sqrt(sos)
	endif  

end subroutine den_ein_2_pre_sos_vandeW

subroutine CVarToPVar_vandeW(nuvar,CVar, PVar , iEOS)
	implicit none
	real*8, dimension(nuvar),intent(in)::CVar
	real*8, dimension(nuvar),intent(inout)::PVar
	real*8  :: gamma4(4),d,p,e,ein,ftau,gg,tau 
	integer :: iEOS,n,nuvar
	if(iEOS < 0 ) stop 'iEOS is Error'

	Call ieos_to_gamma_vandeW(ieos,gamma4)

	n = size(CVar)
	d = CVar(1)
	E = CVar(n)   
	tau = 1/d 
	ein = (E-0.5*dot_product(CVar(2:n-1),CVar(2:n-1))/d) / d

	ftau = gamma4(3)/gamma4(4)*(gamma4(1)/tau/(tau-gamma4(2)))-gamma4(1)/tau/tau 
	gg   = gamma4(3)/gamma4(4)*(tau/(tau-gamma4(2)))
	p = gg * ein * d + ftau 

	PVar(1) = d
	PVar(2:n-1) = CVar(2:n-1)/d
	PVar(n) = p

endsubroutine CVarToPVar_vandeW


subroutine PVarToCVar_vandeW( nUvar, PVar,CVar , iEOS)
	implicit none  
	real*8, dimension(nUvar),intent(inout)::CVar
	real*8, dimension(nUvar),intent(in)::PVar
	real*8  :: gamma4(4),d,p,e,ein,ftau,gg,tau 
	integer :: iEOS,n,nUvar
	if(iEOS < 0 ) stop 'iEOS is Error'   

	Call ieos_to_gamma_vandeW(ieos,gamma4)
	n = size(PVar)
	d = PVar(1)
	p = PVar(n)
	tau = 1/d 
	ftau = gamma4(3)/gamma4(4)*(gamma4(1)/tau/(tau-gamma4(2)))-gamma4(1)/tau/tau 
	gg   = gamma4(3)/gamma4(4)*(tau/(tau-gamma4(2)))
	ein = ( p - ftau )/gg*tau 

	E = 0.5*d*dot_product(PVar(2:n-1),PVar(2:n-1)) + ein*d
	CVar(1) = d
	CVar(2:n-1) = PVar(2:n-1)*d
	CVar(n) = E

endsubroutine PVarToCVar_vandeW

subroutine ieos_to_gamma_vandeW(ieos,gamma4)
	implicit none 
	integer,intent(in)   :: ieos 
	real*8,intent(out)   :: gamma4(4)

	gamma4 = 0.d0 
	!> alpha,beta,R,Cv
	Select case(ieos)
	case(22423)
		gamma4(1) = 0.14d0
		gamma4(2) = 3.258d-5 
		gamma4(3) = 286.9d0
		gamma4(4) = 577.8d0
	case default
		stop 'vandeW EOS is error.'
	End Select

end subroutine ieos_to_gamma_vandeW


subroutine PVarToOther_MieGruSolid( nUvar, PVar,sos,ein, iEOS)
	implicit none  
	real*8, intent(out):: sos,ein
	real*8, dimension(nUvar),intent(in)::PVar
	real*8  :: gamma4(4),d,p,ftau,tau,gg,rho_0,a_0,gma_0,Sm,tau_0
	integer :: iEOS,nUvar,n
	if(iEOS < 0 ) stop 'iEOS is Error'   

	Call ieos_to_gamma_MieGruSolid(ieos,gamma4)

	n = size(PVar)
	d = PVar(1)
	tau = 1/d
	p = PVar(n)

	rho_0 = gamma4(1)
	a_0   = gamma4(2)
	gma_0 = gamma4(3)
	Sm    = gamma4(4)
	tau_0 = 1.d0/rho_0 
	ftau = a_0*a_0*(tau_0-tau)*(tau_0-0.5d0*gma_0*(tau_0-tau)) / &
			tau_0/(tau_0-Sm*(tau_0-tau))**2
	ein = ( p - ftau )/gma_0*tau_0 

	sos = gma_0*p*tau*tau/tau_0 + a_0*a_0*(tau_0/tau-0.5d0*gma_0*(tau_0/tau-1.d0) + &
			(tau_0/tau-1.d0)*(1.d0-0.5d0*gma_0))/(tau_0/tau-sm*(tau_0/tau-1.d0))**2 - &
			2.d0*a_0*a_0*(1.d0-sm)*((tau_0/tau-1.d0)*(tau_0/tau-0.5d0*gma_0*(tau_0/tau-1.d0)))/&
			(tau_0/tau-sm*(tau_0/tau-1.d0))**3

	if ( sos <= 0.d0 ) then 
		print*, sos ,ein
		stop 'sos is less than 0bb'
	else 
		sos = sqrt(sos)
	endif  

endsubroutine PVarToOther_MieGruSolid


subroutine den_ein_2_pre_sos_MieGruSolid(iEOS,d,ein,p,sos)
	implicit none
	integer :: iEOS 
	real*8,intent(in)  :: d,ein
	real*8,intent(out) :: p,sos
	real*8 :: gamma4(4),ftau,gg,tau,rho_0,a_0,gma_0,Sm,tau_0

	Call ieos_to_gamma_MieGruSolid(ieos,gamma4)

	tau = 1/d 
	rho_0 = gamma4(1)
	a_0   = gamma4(2)
	gma_0 = gamma4(3)
	Sm    = gamma4(4)
	tau_0 = 1.d0/rho_0 
	ftau = a_0*a_0*(tau_0-tau)*(tau_0-0.5d0*gma_0*(tau_0-tau)) / &
			tau_0/(tau_0-Sm*(tau_0-tau))**2

	p = gma_0 * ein / tau_0 + ftau 
	sos = gma_0*p*tau*tau/tau_0 + a_0*a_0*(tau_0/tau-0.5d0*gma_0*(tau_0/tau-1.d0) + &
			(tau_0/tau-1.d0)*(1.d0-0.5d0*gma_0))/(tau_0/tau-sm*(tau_0/tau-1.d0))**2 - &
			2.d0*a_0*a_0*(1.d0-sm)*((tau_0/tau-1.d0)*(tau_0/tau-0.5d0*gma_0*(tau_0/tau-1.d0)))/&
			(tau_0/tau-sm*(tau_0/tau-1.d0))**3

	if ( sos <= 0.d0 ) then 
		print*, sos ,ein
		stop 'sos is less than 0bb'
	else 
		sos = sqrt(sos)
	endif  

end subroutine den_ein_2_pre_sos_MieGruSolid

subroutine CVarToPVar_MieGruSolid(nuvar,CVar, PVar , iEOS)
	implicit none
	real*8, dimension(nuvar),intent(in)::CVar
	real*8, dimension(nuvar),intent(inout)::PVar
	real*8  :: gamma4(4),d,p,e,ein,ftau,gg,tau,rho_0,a_0,gma_0,Sm,tau_0
	integer :: iEOS,n,nuvar
	if(iEOS < 0 ) stop 'iEOS is Error'

	Call ieos_to_gamma_MieGruSolid(ieos,gamma4)

	n = size(CVar)
	d = CVar(1)
	E = CVar(n)   
	tau = 1/d 
	ein = (E-0.5*dot_product(CVar(2:n-1),CVar(2:n-1))/d) / d

	rho_0 = gamma4(1)
	a_0   = gamma4(2)
	gma_0 = gamma4(3)
	Sm    = gamma4(4)
	tau_0 = 1.d0/rho_0 
	ftau = a_0*a_0*(tau_0-tau)*(tau_0-0.5d0*gma_0*(tau_0-tau)) / &
			tau_0/(tau_0-Sm*(tau_0-tau))**2

	p = gma_0 * ein /tau_0 + ftau 

	PVar(1) = d
	PVar(2:n-1) = CVar(2:n-1)/d
	PVar(n) = p

endsubroutine CVarToPVar_MieGruSolid


subroutine PVarToCVar_MieGruSolid( nUvar, PVar,CVar , iEOS)
	implicit none  
	real*8, dimension(nUvar),intent(inout)::CVar
	real*8, dimension(nUvar),intent(in)::PVar
	real*8  :: gamma4(4),d,p,e,ein,ftau,gg,tau,rho_0,a_0,gma_0,Sm,tau_0
	integer :: iEOS,n,nUvar
	if(iEOS < 0 ) stop 'iEOS is Error'   

	Call ieos_to_gamma_MieGruSolid(ieos,gamma4)
	n = size(PVar)
	d = PVar(1)
	p = PVar(n)
	tau = 1/d 

	rho_0 = gamma4(1)
	a_0   = gamma4(2)
	gma_0 = gamma4(3)
	Sm    = gamma4(4)
	tau_0 = 1.d0/rho_0 
	ftau = a_0*a_0*(tau_0-tau)*(tau_0-0.5d0*gma_0*(tau_0-tau)) / &
			tau_0/(tau_0-Sm*(tau_0-tau))**2
	ein = ( p - ftau )/gma_0*tau_0 

	E = 0.5*d*dot_product(PVar(2:n-1),PVar(2:n-1)) + ein*d
	CVar(1) = d
	CVar(2:n-1) = PVar(2:n-1)*d
	CVar(n) = E

endsubroutine PVarToCVar_MieGruSolid

subroutine ieos_to_gamma_MieGruSolid(ieos,gamma4)
	implicit none 
	integer,intent(in)   :: ieos 
	real*8,intent(out)   :: gamma4(4)

	gamma4 = 0.d0 
	!> rho_0,a_0,Gamma_0,Sm
	Select case(ieos)
	case(14719)
		gamma4(1) = 2785.d0 
		gamma4(2) = 5328.d0  
		gamma4(3) = 2.d0 
		gamma4(4) = 1.338d0
	case default
		stop 'MieGruSolid EOS is error.'
	End Select

end subroutine ieos_to_gamma_MieGruSolid


End Module mdu_cc
