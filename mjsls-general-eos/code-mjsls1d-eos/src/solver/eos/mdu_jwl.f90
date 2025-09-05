
Module mdu_jwl   
	implicit none

contains   

!> JWL 
subroutine PVarToOther_JWL( nUvar, PVar,sos,ein, iEOS)
	implicit none  
	real*8, intent(out):: sos,ein
	real*8, dimension(nUvar),intent(in)::PVar
	real*8  :: gamma6(6),d,p,ftau
	integer :: iEOS,nUvar,n
	if(iEOS < 0 ) stop 'iEOS is Error'   

	Call ieos_to_gamma_JWL(ieos,gamma6)
	n = size(PVar)
	d = PVar(1)
	p = PVar(n)
	ftau = gamma6(3)*(1.d0-gamma6(1)/gamma6(5)*d/gamma6(2)) * &
				exp(-gamma6(5)*gamma6(2)/d) + &
				gamma6(4)*(1.d0-gamma6(1)/gamma6(6)*d/gamma6(2)) * &
				exp(-gamma6(6)*gamma6(2)/d) 
	ein  = (p - ftau)/gamma6(1)/d
	sos = (gamma6(1)+1.d0)*p/d - &
			gamma6(3)/d*(gamma6(1)+1.d0-gamma6(5)*gamma6(2)/d) &
			* exp(-gamma6(5)*gamma6(2)/d) - & 
			gamma6(4)/d*(gamma6(1)+1.d0-gamma6(6)*gamma6(2)/d) &
			* exp(-gamma6(6)*gamma6(2)/d)	

	if ( sos <= 0.d0 ) then 
		print*, sos ,ein
		stop 'sos is less than 0bb'
	else 
		sos = sqrt(sos)
	endif  

endsubroutine PVarToOther_JWL	


subroutine den_ein_2_pre_sos_JWL(iEOS,d,ein,p,sos)
	implicit none
	integer :: iEOS 
	real*8,intent(in)  :: d,ein
	real*8,intent(out) :: p,sos
	real*8 :: gamma6(6),ftau 

	Call ieos_to_gamma_JWL(ieos,gamma6)
	ftau = gamma6(3)*(1.d0-gamma6(1)/gamma6(5)*d/gamma6(2)) * &
				exp(-gamma6(5)*gamma6(2)/d) + &
				gamma6(4)*(1.d0-gamma6(1)/gamma6(6)*d/gamma6(2)) * &
				exp(-gamma6(6)*gamma6(2)/d) 
	p = ftau + gamma6(1) * d * ein 
	sos = (gamma6(1)+1.d0)*p/d - &
			gamma6(3)/d*(gamma6(1)+1.d0-gamma6(5)*gamma6(2)/d) &
			* exp(-gamma6(5)*gamma6(2)/d) - & 
			gamma6(4)/d*(gamma6(1)+1.d0-gamma6(6)*gamma6(2)/d) &
			* exp(-gamma6(6)*gamma6(2)/d)		
			
	if ( sos <= 0.d0 ) then 
		print*, sos ,ein
		stop 'sos is less than 0bb'
	else 
		sos = sqrt(sos)
	endif  

end subroutine den_ein_2_pre_sos_JWL


subroutine CVarToPVar_JWL(nuvar,CVar, PVar , iEOS)
	implicit none
	real*8, dimension(nuvar),intent(in)::CVar
	real*8, dimension(nuvar),intent(inout)::PVar
	real*8  :: gamma6(6),d,p,e,ein,ftau
	integer :: iEOS,n,nuvar
	if(iEOS < 0 ) stop 'iEOS is Error'

	Call ieos_to_gamma_JWL(ieos,gamma6)

	n = size(CVar)
	d = CVar(1)
	E = CVar(n)   
	ein = (E-0.5*dot_product(CVar(2:n-1),CVar(2:n-1))/d) / d
	ftau = gamma6(3)*(1.d0-gamma6(1)/gamma6(5)*d/gamma6(2)) * &
				exp(-gamma6(5)*gamma6(2)/d) + &
				gamma6(4)*(1.d0-gamma6(1)/gamma6(6)*d/gamma6(2)) * &
				exp(-gamma6(6)*gamma6(2)/d) 
	p = ftau + gamma6(1) * d * ein 
	PVar(1) = d
	PVar(2:n-1) = CVar(2:n-1)/d
	PVar(n) = p

endsubroutine CVarToPVar_JWL


subroutine PVarToCVar_JWL( nUvar, PVar,CVar , iEOS)
	implicit none  
	real*8, dimension(nUvar),intent(inout)::CVar
	real*8, dimension(nUvar),intent(in)::PVar
	real*8  :: gamma6(6),d,p,e,ein,ftau 
	integer :: iEOS,n,nUvar
	if(iEOS < 0 ) stop 'iEOS is Error'   

	Call ieos_to_gamma_JWL(ieos,gamma6)
	n = size(PVar)
	d = PVar(1)
	p = PVar(n)
	ftau = gamma6(3)*(1.d0-gamma6(1)/gamma6(5)*d/gamma6(2)) * &
				exp(-gamma6(5)*gamma6(2)/d) + &
				gamma6(4)*(1.d0-gamma6(1)/gamma6(6)*d/gamma6(2)) * &
				exp(-gamma6(6)*gamma6(2)/d) 
	ein  = (p - ftau)/gamma6(1)/d
	E = 0.5*d*dot_product(PVar(2:n-1),PVar(2:n-1)) + ein*d
	CVar(1) = d
	CVar(2:n-1) = PVar(2:n-1)*d
	CVar(n) = E
endsubroutine PVarToCVar_JWL

subroutine ieos_to_gamma_JWL(ieos,gamma6)
	implicit none 
	integer,intent(in)   :: ieos 
	real*8,intent(out)   :: gamma6(6)

	!> gamma,rho,A,B,R1,R2
	gamma6 = 0.d0 
	Select case(ieos)
	case(17,95022)
		gamma6(1) = 0.8938d0
		gamma6(2) = 1905.d0 
		gamma6(3) = 6.321d13
		gamma6(4) = -4.472d9
		gamma6(5) = 11.3d0 
		gamma6(6) = 1.13d0 
	case(9502,95021)
		gamma6(1) = 0.5d0
		gamma6(2) = 1895.d0 
		gamma6(3) = 1.362d12
		gamma6(4) = 7.199d10
		gamma6(5) = 6.2d0 
		gamma6(6) = 2.2d0 		
	case(3712)
		gamma6(1) = 0.3d0 
		gamma6(2) = 1.63d-3
		gamma6(3) = 3.712d5
		gamma6(4) = 3.23d3
		gamma6(5) = 4.15d0 
		gamma6(6) = 0.95d0 			
	case default
		stop 'JWL EOS is error.'
	End Select
end subroutine ieos_to_gamma_JWL


End Module mdu_jwl
