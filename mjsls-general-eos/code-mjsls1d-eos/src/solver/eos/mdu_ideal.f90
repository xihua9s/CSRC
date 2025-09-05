
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
	case(166,513)
		gamma = 5.d0/3.d0
	case(715)
		gamma = 7.d0 / 5.d0
	case(17,95022)				
		gamma = 0.8938d0	
	case(9502,95021)			
		gamma = 0.5d0		
	case(33)			
		gamma = 1.19d0		
	case(1920)  !> Stiffened gas EOS
		gamma = 4.4d0 	
	case(22423)
		gamma = 0.14d0		
	case(14719) 
		gamma = 2.d0 	
	case(3712)
		gamma = 0.3d0 			
	case default		
		print*, 'iEos = ', ieos 
		stop 'ieos is error, line-gamma'
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


subroutine ieos_to_gamma_stiff(ieos,gamma3)   
	implicit none
	real*8,intent(out) :: gamma3(3)
	integer,intent(in) :: iEOS

	select case(iEos)
	case(1920)  !> Stiffened gas EOS
		gamma3(1) = 	4.4d0
		gamma3(2) = 	0.d0    !> q_infty
		gamma3(3) = 	6.d8	!> pi_infty 
	case(7153309)  !> Stiffened gas EOS
		gamma3(1) = 	7.15d0
		gamma3(2) = 	0.d0    !> q_infty
		gamma3(3) = 	3.309d2	!> pi_infty 
	case default		
		stop 'ieos is error, line11'
	end select

end subroutine ieos_to_gamma_stiff


subroutine PvarToCvar_Stiff( nUvar, PVar,CVar, iEOS)
	implicit none  
	real*8, dimension(nUvar),intent(in)  :: PVar
	real*8, dimension(nUvar),intent(out) :: CVar
	real*8  :: gamma3(3),d,p
	integer :: iEOS,nUvar,n,ein 
	if(iEOS < 0 ) stop 'iEOS is Error'   

	Call ieos_to_gamma_stiff(ieos,gamma3)
	n = size(PVar)
	d = PVar(1)
	p = PVar(n)
	ein = ((p+gamma3(1)*gamma3(3))/(gamma3(1)-1.d0)+d*gamma3(2))/d
	Cvar(1) = d 
	Cvar(2:n-1) = PVar(2:n-1)*d
	Cvar(n) = d*ein + 0.5*d*dot_product(PVar(2:n-1),PVar(2:n-1))

endsubroutine PvarToCvar_Stiff


subroutine PvarToOther_Stiff( nUvar, PVar,sos,ein, iEOS)
	implicit none  
	real*8, intent(out):: sos,ein
	real*8, dimension(nUvar),intent(in)::PVar
	real*8  :: gamma3(3),d,p,ftau
	integer :: iEOS,nUvar,n
	if(iEOS < 0 ) stop 'iEOS is Error'   

	Call ieos_to_gamma_stiff(ieos,gamma3)
	n = size(PVar)
	d = PVar(1)
	p = PVar(n)
	ein = ((p+gamma3(1)*gamma3(3))/(gamma3(1)-1.d0)+d*gamma3(2))/d
	sos = gamma3(1)*(p+gamma3(3))/d

	if ( sos <= 0.d0 ) then 
		print*, sos ,ein
		stop 'sos is less than 0bb'
	else 
		sos = sqrt(sos)
	endif  

endsubroutine PvarToOther_Stiff


subroutine den_ein_2_pre_sos_stiff(iEOS,d,ein,p,sos)
	implicit none
	integer :: iEOS 
	real*8,intent(in)  :: d,ein
	real*8,intent(out) :: p,sos
	real*8 :: gamma3(3) 

	Call ieos_to_gamma_stiff(ieos,gamma3)

	p = (gamma3(1)-1.d0)*(d*ein-d*gamma3(2))-gamma3(1)*gamma3(3)
	sos = gamma3(1)*(p+gamma3(3))/d	
			
	if ( sos <= 0.d0 ) then 
		print*, sos ,ein
		stop 'sos is less than 0bb'
	else 
		sos = sqrt(sos)
	endif  

end subroutine den_ein_2_pre_sos_stiff


End Module mdu_ideal
