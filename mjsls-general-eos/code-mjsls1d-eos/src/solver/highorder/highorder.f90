!> Last Modified : Wed May  4 18:30:54 2022

   !real*8  :: rho3,sos3,vel3,pre3
   !real*8  :: rho4,sos4,vel4,pre4
   !real*8  :: rho_left_slope,rho_rght_slope
   !real*8  :: sos_left_slope,sos_rght_slope
   !real*8  :: vel_left_slope,vel_rght_slope
   !real*8  :: pre_left_slope,pre_rght_slope

   !> 二阶精度
      !do ix = 3 , ini_cmt%nx-1
      !
      !   x1 = mesh_phys%vertex(ix-2)
      !   x2 = mesh_phys%vertex(ix-1)
      !   x3 = mesh_phys%vertex(ix-0)
      !   x4 = mesh_phys%vertex(ix+1)
      !   x5 = mesh_phys%vertex(ix+2)
      !            
      !   rho1 = mesh_phys%con_phys(ix-2,1)
      !   sos1 = mesh_phys%sos(ix-2)
      !   vel1 = mesh_phys%con_phys(ix-2,2) / mesh_phys%con_phys(ix-2,1)
      !   pre1 = mesh_phys%pre(ix-2)
      !
      !   rho2 = mesh_phys%con_phys(ix-1,1)
      !   sos2 = mesh_phys%sos(ix-1)
      !   vel2 = mesh_phys%con_phys(ix-1,2) / mesh_phys%con_phys(ix-1,1)
      !   pre2 = mesh_phys%pre(ix-1)
      !
      !   rho3 = mesh_phys%con_phys(ix-0,1)
      !   sos3 = mesh_phys%sos(ix-0)
      !   vel3 = mesh_phys%con_phys(ix-0,2) / mesh_phys%con_phys(ix-0,1)
      !   pre3 = mesh_phys%pre(ix-0)
      !
      !   rho4 = mesh_phys%con_phys(ix+1,1)
      !   sos4 = mesh_phys%sos(ix+1)
      !   vel4 = mesh_phys%con_phys(ix+1,2) / mesh_phys%con_phys(ix+1,1)
      !   pre4 = mesh_phys%pre(ix+1)
      !
      !   Call slope_1D_cell_centered(x1,x2,x3,x4,rho1,rho2,rho3,rho_left_slope)
      !   Call slope_1D_cell_centered(x2,x3,x4,x5,rho2,rho3,rho4,rho_rght_slope)
      !
      !   Call slope_1D_cell_centered(x1,x2,x3,x4,sos1,sos2,sos3,sos_left_slope)
      !   Call slope_1D_cell_centered(x2,x3,x4,x5,sos2,sos3,sos4,sos_rght_slope)
      !
      !   Call slope_1D_cell_centered(x1,x2,x3,x4,vel1,vel2,vel3,vel_left_slope)
      !   Call slope_1D_cell_centered(x2,x3,x4,x5,vel2,vel3,vel4,vel_rght_slope)
      !
      !   Call slope_1D_cell_centered(x1,x2,x3,x4,pre1,pre2,pre3,pre_left_slope)
      !   Call slope_1D_cell_centered(x2,x3,x4,x5,pre2,pre3,pre4,pre_rght_slope)
      !
      !   a1 = (rho2 + (x3-x2)/2.d0 * rho_left_slope) * (sos2 + (x3-x2)/2.d0 * sos_left_slope) 
      !   r1 = (rho2 + (x3-x2)/2.d0 * rho_left_slope)
      !   u1 = (vel2 + (x3-x2)/2.d0 * vel_left_slope)
      !   p1 = (pre2 + (x3-x2)/2.d0 * pre_left_slope)
      !   
      !   ! a1 = mesh_phys%con_phys(ix-1,1) * mesh_phys%sos(ix-1)
      !   ! u1 = mesh_phys%con_phys(ix-1,2)/mesh_phys%con_phys(ix-1,1)
      !   ! p1 = mesh_phys%pre(ix-1)
      !   
      !   a2 = (rho3 + (x3-x4)/2.d0 * rho_rght_slope) * (sos3 + (x3-x4)/2.d0 * sos_rght_slope)
      !   r2 = (rho3 + (x3-x4)/2.d0 * rho_rght_slope)
      !   u2 = (vel3 + (x3-x4)/2.d0 * vel_rght_slope)
      !   p2 = (pre3 + (x3-x4)/2.d0 * pre_rght_slope)
      !   
      !   ! a2 = mesh_phys%con_phys(ix,1) * mesh_phys%sos(ix)
      !   ! u2 = mesh_phys%con_phys(ix,2)/mesh_phys%con_phys(ix,1)
      !   ! p2 = mesh_phys%pre(ix)
      !   
      !   !> 弱波近似
      !   !mesh_phys%vel_p(ix) = (a1*u1+a2*u2)/(a1+a2) + (p1-p2)/(a1+a2)
      !   !mesh_phys%pstar(ix) = (p1/a1+p2/a2+(u1-u2)) / (1.d0/a1 + 1.d0/a2)
      !   
      !   !> 精确RP解
      !   call ieos_to_gamma(mesh_phys%eos(ix),gamma)
      !   Call exact_riemann_problem_pstar_ustar(r1,u1,p1,r2,u2,p2,gamma,pstar,ustar)
      !   mesh_phys%vel_p(ix) = ustar
      !   mesh_phys%pstar(ix) = pstar        
      !   
      !enddo   


subroutine slope_1D_cell_centered(x1,x2,x3,x4,rho1,rho2,rho3,rho_slope)
   implicit none
   real*8,intent(in)  :: x1,x2,x3,x4,rho1,rho2,rho3
   real*8,intent(out) :: rho_slope

   real*8 :: slope1,slope2,slope3,s,s1,s2,s3
   
   rho_slope = 0.d0
   
   slope1 = (rho2 - rho1) / ( (x3+x2)/2.d0 - (x2+x1)/2.d0 )
   slope2 = (rho3 - rho2) / ( (x4+x3)/2.d0 - (x3+x2)/2.d0 )
   slope3 = (rho3 - rho1) / ( (x4+x3)/2.d0 - (x2+x1)/2.d0 )
   
   
   if( slope1 > 0.d0 .and. slope2 > 0.d0 ) then
      rho_slope = min(slope1,slope2)
   elseif( slope1 < 0.d0 .and. slope2 < 0.d0 ) then
      rho_slope = max(slope1,slope2)
   else
      rho_slope = 0.d0
   endif
   
   rho_slope = 0.d0
   
   !> 选取绝对值最小的那个
   !s1 = abs(slope1)
   !s2 = abs(slope2)
   !s3 = abs(slope3)

   !s = min(s1,min(s2,s3))
   !if(abs(s-s1)<1.d-14) then
   !   rho_slope = slope1
   !elseif(abs(s-s2)<1.d-14) then
   !   rho_slope = slope2
   !elseif(abs(s-s3)<1.d-14) then
   !   rho_slope = slope3
   !else
   !   print*, slope1,slope2,slope3
   !   stop 'slope_1D_cell_centered is Error, pls Call Xihua'
   !endif
   
   !rho_slope = 0.d0
   
   
   
   
end subroutine slope_1D_cell_centered
