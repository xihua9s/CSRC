
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



End Module mdu_cc
