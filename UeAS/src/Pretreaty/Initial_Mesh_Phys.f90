


Subroutine Initial_Mesh_Phys(Ini_CDT,MeshCY2D,PhysCY2D)

   use COMMON_TYPE_Condition
   use COMMON_TYPE_Mesh
   use COMMON_TYPE_Phys
   use mdu_sod_2D_cat_cyl
   implicit none
   
   Type(COMMON_TYPE_CDT),intent(inout) :: Ini_CDT
   Type(Mesh_Cyl_2D),intent(inout)  :: MeshCY2D
   Type(Phys_Cyl_2D),intent(inout)  :: PhysCY2D

   Integer(KNDI)  :: ix,iy
   print*, '---> Initial_Mesh_Phys'
   
   Select Case(Ini_CDT%Test_Case)
   
   Case('Sod')

      Call Initial_Mesh_Phys_Sod(Ini_CDT,MeshCY2D,PhysCY2D)
            
   Case default
   
      Stop 'Ini_CDT%Test_Case is Error, pls Call Xihua.'
   
   End Select
   
   Call Cal_Geometry(Ini_CDT,MeshCY2D)
   
   Do ix = 1 , Ini_CDT%nx
   Do iy = 1 , Ini_CDT%ny
      PhysCY2D%mass(ix,iy) = PhysCY2D%den(ix,iy) * MeshCY2D%volm(ix,iy)
   Enddo
   Enddo

   print*, '<--- Initial_Mesh_Phys'
   
End Subroutine Initial_Mesh_Phys



Subroutine Register_Mesh_Phys(Ini_CDT,MeshCY2D,PhysCY2D)

   use COMMON_TYPE_Condition
   use COMMON_TYPE_Mesh
   use COMMON_TYPE_Phys
   use COMMON_XXX_Constant
   implicit none
   
   Type(COMMON_TYPE_CDT),intent(in) :: Ini_CDT
   Type(Mesh_Cyl_2D),intent(out)    :: MeshCY2D
   Type(Phys_Cyl_2D),intent(out)    :: PhysCY2D
   
   integer(KNDI)  :: nx,ny
   
   print*, '---> Register_Mesh_Phys(Ini_CDT,MeshCY2D,PhysCY2D)'   
   
   nx = Ini_CDT%nx
   ny = Ini_CDT%ny
   
   allocate(PhysCY2D%den(nx,ny)) ; PhysCY2D%den = lageps
   allocate(PhysCY2D%pre(nx,ny)) ; PhysCY2D%pre = lageps
   allocate(PhysCY2D%tau(nx,ny)) ; PhysCY2D%tau = lageps
   allocate(PhysCY2D%enin(nx,ny)) ; PhysCY2D%enin = lageps
   allocate(PhysCY2D%ener(nx,ny)) ; PhysCY2D%ener = lageps
   allocate(PhysCY2D%mass(nx,ny)) ; PhysCY2D%mass = lageps
   allocate(PhysCY2D%sound(nx,ny)) ; PhysCY2D%sound = lageps 
   allocate(PhysCY2D%vel_x(nx,ny)) ; PhysCY2D%vel_x = lageps
   allocate(PhysCY2D%vel_y(nx,ny)) ; PhysCY2D%vel_y = lageps
   allocate(PhysCY2D%EOS(nx,ny))   ; PhysCY2D%EOS = ILageps
   
   allocate(PhysCY2D%vstar_x(0:nx,0:ny)) ; PhysCY2D%vstar_x = lageps 
   allocate(PhysCY2D%vstar_y(0:nx,0:ny)) ; PhysCY2D%vstar_y = lageps
   allocate(PhysCY2D%PstarC1(0:nx,0:ny)) ; PhysCY2D%PstarC1 = lageps
   allocate(PhysCY2D%Pstar8(nx,ny,8)) ; PhysCY2D%Pstar8 = lageps
   allocate(PhysCY2D%Pstar4(nx,ny,4)) ; PhysCY2D%Pstar4 = lageps
   allocate(PhysCY2D%PstarC4(nx,ny,4)) ; PhysCY2D%PstarC4 = lageps
   allocate(PhysCY2D%PstarC8(nx,ny,8)) ; PhysCY2D%PstarC4 = lageps

   allocate(MeshCY2D%area(nx,ny)) ; MeshCY2D%area = lageps
   allocate(MeshCY2D%volm(nx,ny)) ; MeshCY2D%volm = lageps
   allocate(MeshCY2D%vertex_x(0:nx,0:ny)) ; MeshCY2D%vertex_x = lageps
   allocate(MeshCY2D%vertex_y(0:nx,0:ny)) ; MeshCY2D%vertex_y = lageps
   
   allocate(MeshCY2D%node_type_char(0:nx,0:ny))       ; MeshCY2D%node_type_char      = 'inner'
   allocate(MeshCY2D%node_type_direction(0:nx,0:ny,2)); MeshCY2D%node_type_direction = 0.d0
   allocate(MeshCY2D%node_type_velocity(0:nx,0:ny,2)) ; MeshCY2D%node_type_velocity  = 0.d0 
   allocate(MeshCY2D%node_type_pressure(0:nx,0:ny))   ; MeshCY2D%node_type_pressure  = 0.d0 
   
   print*, '<--- Register_Mesh_Phys(Ini_CDT,MeshCY2D,PhysCY2D)'   
   
End Subroutine Register_Mesh_Phys



