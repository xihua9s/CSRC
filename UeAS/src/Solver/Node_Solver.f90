
Subroutine Node_Solver(Ini_CDT,MeshCY2D,PhysCY2D)
   use COMMON_TYPE_Condition
   use COMMON_TYPE_Mesh
   use COMMON_TYPE_Phys
   use COMMON_XXX_Constant
   implicit none

   Type(COMMON_TYPE_CDT),intent(inout) :: Ini_CDT
   Type(Mesh_Cyl_2D),intent(inout) :: MeshCY2D
   Type(Phys_Cyl_2D),intent(inout) :: PhysCY2D

   Select Case(Ini_CDT%Node_Solver)
      
   Case('UeAS')

      Select Case(Ini_CDT%Cat_Cyl)
      
      Case('Cat')

         Call Node_Solver_Cat_UeAS(Ini_CDT,MeshCY2D,PhysCY2D)
         
      End Select  

   Case default

   stop 'Ini_CDT%Node_Solver is Error, pls Call Xihua.'

   End Select   

End Subroutine Node_Solver



