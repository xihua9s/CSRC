


Subroutine Solve_EulerEqn(Ini_CDT,MeshCY2D,PhysCY2D)
   use COMMON_TYPE_Condition
   use COMMON_TYPE_Mesh
   use COMMON_TYPE_Phys
   use COMMON_XXX_Constant
   implicit none

   Type(COMMON_TYPE_CDT),intent(inout) :: Ini_CDT
   Type(Mesh_Cyl_2D),intent(inout) :: MeshCY2D
   Type(Phys_Cyl_2D),intent(inout) :: PhysCY2D

   Call Output_initial_Mesh_Phys(Ini_CDT,MeshCY2D,PhysCY2D)   

   Do While(Ini_CDT%time < Ini_CDT%Tstop)

      Call Output_Mesh_Phys(Ini_CDT,MeshCY2D,PhysCY2D)

      Call Get_Dt(Ini_CDT,MeshCY2D,PhysCY2D)

      Call Node_Solver(Ini_CDT,MeshCY2D,PhysCY2D)

      Call Update_Mesh_Phys(Ini_CDT,MeshCY2D,PhysCY2D)
           
      if(Ini_CDT%it > Ini_CDT%It_num) stop 'Ini_CDT%it > Ini_CDT%It_num'

   Enddo

   Call Output_final_Mesh_Phys(Ini_CDT,MeshCY2D,PhysCY2D)   

End Subroutine Solve_EulerEqn
   