
Subroutine Test_Benchmarck_Case(Ini_CDT)

   use COMMON_TYPE_Condition
   use COMMON_TYPE_Mesh
   use COMMON_TYPE_Phys
   implicit none
   
   Type(COMMON_TYPE_CDT) :: Ini_CDT
   Type(Mesh_Cyl_2D) :: MeshCY2D
   Type(Phys_Cyl_2D) :: PhysCY2D

   print*, '---> Test_Benchmarck_Case'
   
   Call Register_Mesh_Phys(Ini_CDT,MeshCY2D,PhysCY2D)
 
   Call Initial_Mesh_Phys(Ini_CDT,MeshCY2D,PhysCY2D)
   
   Call Solve_EulerEqn(Ini_CDT,MeshCY2D,PhysCY2D)
  
   Call Output_Mesh_Phys(Ini_CDT,MeshCY2D,PhysCY2D)
  
   print*, '<--- Test_Benchmarck_Case'
   
End Subroutine Test_Benchmarck_Case  

