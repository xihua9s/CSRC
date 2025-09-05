
Module mdu_ideal   
	implicit none

contains   

subroutine ieos_to_gamma(ieos,gamma)   
	implicit none
	real*8,intent(out) :: gamma
	integer,intent(in) :: iEOS

	select case(iEos)
	case(140)
		gamma = 1.4d0
	case(166)
		gamma = 5.d0/3.d0
	case(715)
		gamma = 7.d0 / 5.d0
	case(513)
		gamma = 5.d0 / 3.d0     
	case(17)				
		gamma = 0.8938d0	
	case(9502)			
		gamma = 0.5d0		
	case(95021)			
		gamma = 0.5d0		
	case(95022)				
		gamma = 0.8938d0	
	case(33)			
		gamma = 1.19d0		
	case default		
		stop 'ieos is error, line11'
	end select

end subroutine ieos_to_gamma

subroutine CVarToPVar(nuvar,CVar, PVar , gamma , iEOS)
	implicit none
	real*8, dimension(nuvar),intent(in)::CVar
	real*8, dimension(nuvar),intent(inout)::PVar
	real*8  :: gamma,d,e,p
	integer :: iEOS,n,nuvar
	if(iEOS < 0 ) stop 'iEOS is Error'

	n = size(CVar)
	d = CVar(1)
	E = CVar(n)   
	p = (E-0.5*dot_product(CVar(2:n-1),CVar(2:n-1))/d)*(gamma-1.d0)
	PVar(1) = d
	PVar(2:n-1) = CVar(2:n-1)/d
	PVar(n) = p

endsubroutine CVarToPVar

!   
subroutine PVarToCVar( nUvar, PVar,CVar, Gamma , iEOS)
	implicit none  
	real*8, dimension(nUvar),intent(inout)::CVar
	real*8, dimension(nUvar),intent(in)::PVar
	real*8  :: gamma,d,p,e
	integer :: iEOS,n,nUvar
	if(iEOS < 0 ) stop 'iEOS is Error'   

	n = size(PVar)
	d = PVar(1)
	p = PVar(n)
	E = 0.5*d*dot_product(PVar(2:n-1),PVar(2:n-1))+p/(Gamma-1.d0)
	CVar(1) = d
	CVar(2:n-1) = PVar(2:n-1)*d
	CVar(n) = E
endsubroutine PVarToCVar


subroutine PVarToOther( nUvar, PVar,sos,ein, Gamma , iEOS)
	implicit none  
	real*8, intent(out):: sos,ein
	real*8, dimension(nUvar),intent(in)::PVar
	real*8  :: gamma
	integer :: iEOS,nUvar
	if(iEOS < 0 ) stop 'iEOS is Error'   

	sos = sqrt(gamma*PVar(nUvar)/PVar(1))      
	ein = PVar(nUvar)/PVar(1)/(gamma-1.d0)

endsubroutine PVarToOther


End Module mdu_ideal
