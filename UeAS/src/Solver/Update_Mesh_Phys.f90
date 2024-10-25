
Subroutine Update_Mesh_Phys(Ini_CDT,MeshCY2D,PhysCY2D)

   use COMMON_TYPE_Condition
   use COMMON_TYPE_Mesh
   use COMMON_TYPE_Phys
   implicit none
   
   Type(COMMON_TYPE_CDT),intent(in) :: Ini_CDT
   Type(Mesh_Cyl_2D),intent(inout) :: MeshCY2D
   Type(Phys_Cyl_2D),intent(inout) :: PhysCY2D
   
   Select Case(Ini_CDT%Cat_Cyl)
   
   Case('Cat')
      
      Call Update_Mesh_Phys_Cat(Ini_CDT,MeshCY2D,PhysCY2D)
      
   End Select

End Subroutine Update_Mesh_Phys

Subroutine Update_Mesh_Phys_Cat(Ini_CDT,MeshCY2D,PhysCY2D)

   use COMMON_TYPE_Condition
   use COMMON_TYPE_Mesh
   use COMMON_TYPE_Phys
   implicit none
   
   Type(COMMON_TYPE_CDT),intent(in) :: Ini_CDT
   Type(Mesh_Cyl_2D),intent(inout) :: MeshCY2D
   Type(Phys_Cyl_2D),intent(inout) :: PhysCY2D
   
   Select Case(Ini_CDT%Node_Solver)
   
   Case('UeAS')
   
      Select Case(Ini_CDT%Origin_YesNo)
      
      Case('Yes') 

         Call Update_Mesh_Phys_Cat_Xihua_0_1017(PhysCY2D%den,PhysCY2D%vel_x,PhysCY2D%vel_y,PhysCY2D%pre,    &
            PhysCY2D%sound,PhysCY2D%tau,PhysCY2D%enin,PhysCY2D%ener,PhysCY2D%eos,PhysCY2D%vstar_x,&
            PhysCY2D%vstar_y,PhysCY2D%PstarC1,PhysCY2D%mass,MeshCY2D%vertex_x,MeshCY2D%vertex_y,    &
            MeshCY2D%area,MeshCY2D%volm,Ini_CDT%nx,Ini_CDT%ny,Ini_CDT%dt)       
                           
      Case default
      
         print*, 'Ini_CDT%Origin_YesNo is:', Ini_CDT%Origin_YesNo
         stop 'Error, fdajl39fdja'
      
      End Select
      
   Case default
   
      stop 'Ini_CDT%Node_Solver is Error, pls Call Xihua.'
      
   End Select
      
   Call Cal_Geometry(Ini_CDT,MeshCY2D)        

End Subroutine Update_Mesh_Phys_Cat
