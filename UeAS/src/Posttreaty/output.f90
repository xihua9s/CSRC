
Subroutine Output_initial_Mesh_Phys(Ini_CDT,MeshCY2D,PhysCY2D)   
   use COMMON_TYPE_Condition
   use COMMON_TYPE_Mesh
   use COMMON_TYPE_Phys
   implicit none
   
   Type(COMMON_TYPE_CDT),intent(in) :: Ini_CDT
   Type(Mesh_Cyl_2D) :: MeshCY2D
   Type(Phys_Cyl_2D) :: PhysCY2D
   
   Call output_initial(Ini_CDT%nx,Ini_CDT%ny,PhysCY2D%den,PhysCY2D%vel_x,PhysCY2D%vel_y, & 
         PhysCY2D%pre,PhysCY2D%enin,PhysCY2D%ener,PhysCY2D%mass,MeshCY2D%vertex_x,       &
         MeshCY2D%vertex_y,Ini_CDT%it,Ini_CDT%time)
   
End Subroutine Output_initial_Mesh_Phys

Subroutine Output_final_Mesh_Phys(Ini_CDT,MeshCY2D,PhysCY2D)   
   use COMMON_TYPE_Condition
   use COMMON_TYPE_Mesh
   use COMMON_TYPE_Phys
   implicit none
   
   Type(COMMON_TYPE_CDT),intent(in) :: Ini_CDT
   Type(Mesh_Cyl_2D) :: MeshCY2D
   Type(Phys_Cyl_2D) :: PhysCY2D
   
   Call output_final(Ini_CDT%nx,Ini_CDT%ny,PhysCY2D%den,PhysCY2D%vel_x,PhysCY2D%vel_y, & 
         PhysCY2D%pre,PhysCY2D%enin,PhysCY2D%ener,PhysCY2D%mass,MeshCY2D%vertex_x,       &
         MeshCY2D%vertex_y,Ini_CDT%it,Ini_CDT%time)
   
End Subroutine Output_final_Mesh_Phys

Subroutine Output_Mesh_Phys(Ini_CDT,MeshCY2D,PhysCY2D)   
   use COMMON_TYPE_Condition
   use COMMON_TYPE_Mesh
   use COMMON_TYPE_Phys
   implicit none
   
   Type(COMMON_TYPE_CDT),intent(in) :: Ini_CDT
   Type(Mesh_Cyl_2D) :: MeshCY2D
   Type(Phys_Cyl_2D) :: PhysCY2D
   
   if(mod(Ini_CDT%it,Ini_CDT%It_out)==0) then
      print*, '===> Output_Mesh_Phys, num. is: ', Ini_CDT%it
      Call output(Ini_CDT%nx,Ini_CDT%ny,PhysCY2D%den,PhysCY2D%vel_x,PhysCY2D%vel_y,       & 
            PhysCY2D%pre,PhysCY2D%enin,PhysCY2D%ener,PhysCY2D%mass,MeshCY2D%vertex_x,     &
            MeshCY2D%vertex_y,Ini_CDT%it,Ini_CDT%time)
   endif
   
End Subroutine Output_Mesh_Phys


Subroutine output_initial(nx,ny,den,velx,vely,pre,enin,ener,mass,vertx,verty,it,time)
   use COMMON_XXX_Constant
   implicit none

   integer,intent(in) :: nx,ny,it
   real*8,intent(in)  :: time
   real*8,intent(in)  :: den(nx,ny),velx(nx,ny),vely(nx,ny)
   real*8,intent(in)  :: pre(nx,ny),enin(nx,ny),ener(nx,ny),mass(nx,ny)
   real*8,intent(in)  :: vertx(0:nx,0:ny),verty(0:nx,0:ny)
   
   integer :: i,j
   real*8  :: xmid,ymid,r
   character(99) :: fname

   write(fname,'(i10.10)') it
   open(unit=100,file=FoldOutput//'/phys_initial.plt')
   write(100,*) 'title="contour"'
   write(100,*) 'variables="x","y","den","pre","enin","ener","mass","velx","vely"'
   write(100,*) 'zone  solutiontime=',time, 'i=', nx+1, 'j=' ,ny+1
   write(100,*) 'DATAPACKING=BLOCK,VARLOCATION=([3-9]=CELLCENTERED)'
   do j = 0 , ny 
   do i = 0 , nx 
      WRITE(100,101)  vertx(i,j)
   enddo
   enddo
   do j = 0 , ny
   do i = 0 , nx
      WRITE(100,101)  verty(i,j)
   enddo
   enddo
   do j = 1 , ny
   do i = 1 , nx
      WRITE(100,101)  den(i,j)
   enddo
   enddo   
   do j = 1 , ny
   do i = 1 , nx
      WRITE(100,101)  pre(i,j)
   enddo
   enddo   
   do j = 1 , ny
   do i = 1 , nx
      WRITE(100,101)  enin(i,j)
   enddo
   enddo   
   do j = 1 , ny
   do i = 1 , nx
      WRITE(100,101)  ener(i,j)
   enddo
   enddo   
   do j = 1 , ny
   do i = 1 , nx
      WRITE(100,101)  mass(i,j)
   enddo
   enddo   
   do j = 1 , ny
   do i = 1 , nx
      WRITE(100,101)  velx(i,j)
   enddo
   enddo   
   do j = 1 , ny
   do i = 1 , nx
      WRITE(100,101)  vely(i,j)
   enddo
   enddo   
   
101 format((F14.6))
   close(100)

end subroutine output_initial


Subroutine output_final(nx,ny,den,velx,vely,pre,enin,ener,mass,vertx,verty,it,time)
   use COMMON_XXX_Constant
   implicit none

   integer,intent(in) :: nx,ny,it
   real*8,intent(in)  :: time
   real*8,intent(in)  :: den(nx,ny),velx(nx,ny),vely(nx,ny)
   real*8,intent(in)  :: pre(nx,ny),enin(nx,ny),ener(nx,ny),mass(nx,ny)
   real*8,intent(in)  :: vertx(0:nx,0:ny),verty(0:nx,0:ny)
   
   integer :: i,j
   real*8  :: xmid,ymid,r,r1,r2,r3,r4
   character(99) :: fname

   write(fname,'(i10.10)') it
   open(unit=100,file=FoldOutput//'/phys_final.plt')
   write(100,*) 'title="contour"'
   write(100,*) 'variables="x","y","den","pre","enin","ener","mass","velx","vely"'
   write(100,*) 'zone  solutiontime=',time, 'i=', nx+1, 'j=' ,ny+1
   write(100,*) 'DATAPACKING=BLOCK,VARLOCATION=([3-9]=CELLCENTERED)'
   do j = 0 , ny 
   do i = 0 , nx 
      WRITE(100,101)  vertx(i,j)
   enddo
   enddo
   do j = 0 , ny
   do i = 0 , nx
      WRITE(100,101)  verty(i,j)
   enddo
   enddo
   do j = 1 , ny
   do i = 1 , nx
      WRITE(100,101)  den(i,j)
   enddo
   enddo   
   do j = 1 , ny
   do i = 1 , nx
      WRITE(100,101)  pre(i,j)
   enddo
   enddo   
   do j = 1 , ny
   do i = 1 , nx
      WRITE(100,101)  enin(i,j)
   enddo
   enddo   
   do j = 1 , ny
   do i = 1 , nx
      WRITE(100,101)  ener(i,j)
   enddo
   enddo   
   do j = 1 , ny
   do i = 1 , nx
      WRITE(100,101)  mass(i,j)
   enddo
   enddo   
   do j = 1 , ny
   do i = 1 , nx
      WRITE(100,101)  velx(i,j)
   enddo
   enddo   
   do j = 1 , ny
   do i = 1 , nx
      WRITE(100,101)  vely(i,j)
   enddo
   enddo   
   
101 format((F14.6))
   close(100)
   
   open(unit=100,file=FoldOutput//'/phys_final_R.plt')
   write(100,*) 'VARIABLES = "X","Density","Pressure","Velocity","Internal Energy"'
   write(100,*)  'ZONE I =', nx, ', F=POINT' 
   
   do j = ny/2 , ny/2
   do i = 1 , nx
      r1 = sqrt(vertx(i-1,j)*vertx(i-1,j) + verty(i-1,j)*verty(i-1,j))
      r2 = sqrt(vertx(i-1,j-1)*vertx(i-1,j-1) + verty(i-1,j-1)*verty(i-1,j-1))
      r3 = sqrt(vertx(i,j-1)*vertx(i,j-1) + verty(i,j-1)*verty(i,j-1))
      r4 = sqrt(vertx(i,j)*vertx(i,j) + verty(i,j)*verty(i,j))
      r = (r1+r2+r3+r4)/4.d0
      WRITE(100,102)  r , den(i,j) , pre(i,j) , sqrt(velx(i,j)*velx(i,j)+vely(i,j)*vely(i,j)),enin(i,j)
   enddo
   enddo      
102 format(5(F14.6))
   close(100)   
   
end subroutine output_final

Subroutine output(nx,ny,den,velx,vely,pre,enin,ener,mass,vertx,verty,it,time)
   use COMMON_XXX_Constant
   implicit none

   integer,intent(in) :: nx,ny,it
   real*8,intent(in)  :: time
   real*8,intent(in)  :: den(nx,ny),velx(nx,ny),vely(nx,ny)
   real*8,intent(in)  :: pre(nx,ny),enin(nx,ny),ener(nx,ny),mass(nx,ny)
   real*8,intent(in)  :: vertx(0:nx,0:ny),verty(0:nx,0:ny)
   
   integer :: i,j
   real*8  :: xmid,ymid,r,r1,r2,r3,r4
   character(99) :: fname
   
   write(fname,'(i10.10)') it
   open(unit=100,file=FoldOutput//'/phys_'//trim(adjustl(fname))//'.plt')
   write(100,*) 'title="contour"'
   write(100,*) 'variables="x","y","den","pre","enin","ener","mass","velx","vely"'
   write(100,*) 'zone  solutiontime=',time, 'i=', nx+1, 'j=' ,ny+1
   write(100,*) 'DATAPACKING=BLOCK,VARLOCATION=([3-9]=CELLCENTERED)'
   do j = 0 , ny 
   do i = 0 , nx 
      WRITE(100,101)  vertx(i,j)
   enddo
   enddo
   do j = 0 , ny
   do i = 0 , nx
      WRITE(100,101)  verty(i,j)
   enddo
   enddo
   do j = 1 , ny
   do i = 1 , nx
      WRITE(100,101)  den(i,j)
   enddo
   enddo   
   do j = 1 , ny
   do i = 1 , nx
      WRITE(100,101)  pre(i,j)
   enddo
   enddo   
   do j = 1 , ny
   do i = 1 , nx
      WRITE(100,101)  enin(i,j)
   enddo
   enddo   
   do j = 1 , ny
   do i = 1 , nx
      WRITE(100,101)  ener(i,j)
   enddo
   enddo   
   do j = 1 , ny
   do i = 1 , nx
      WRITE(100,101)  mass(i,j)
   enddo
   enddo   
   do j = 1 , ny
   do i = 1 , nx
      WRITE(100,101)  velx(i,j)
   enddo
   enddo   
   do j = 1 , ny
   do i = 1 , nx
      WRITE(100,101)  vely(i,j)
   enddo
   enddo   
   
101 format((F14.6))
   close(100)

   open(unit=100,file=FoldOutput//'/phys_'//trim(adjustl(fname))//'_R.plt')      
   write(100,*) 'VARIABLES = "X","Density","Pressure","Velocity","Internal Energy"'
   write(100,*)  'ZONE I =', nx, ', F=POINT' 
   
   do j = ny/2 , ny/2
   do i = 1 , nx
      r1 = sqrt(vertx(i-1,j)*vertx(i-1,j) + verty(i-1,j)*verty(i-1,j))
      r2 = sqrt(vertx(i-1,j-1)*vertx(i-1,j-1) + verty(i-1,j-1)*verty(i-1,j-1))
      r3 = sqrt(vertx(i,j-1)*vertx(i,j-1) + verty(i,j-1)*verty(i,j-1))
      r4 = sqrt(vertx(i,j)*vertx(i,j) + verty(i,j)*verty(i,j))
      r = (r1+r2+r3+r4)/4.d0
      WRITE(100,102)  r , den(i,j) , pre(i,j) , sqrt(velx(i,j)*velx(i,j)+vely(i,j)*vely(i,j)),enin(i,j)
   enddo
   enddo      
102 format(5(F14.6))
   close(100)      
   
end subroutine output