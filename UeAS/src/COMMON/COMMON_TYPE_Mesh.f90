


Module COMMON_TYPE_Mesh
   use COMMON_const_kind
   implicit none
   
   Type :: Mesh_Cyl_2D
   
      real(KNDR),allocatable  ::  vertex_x(:,:)
      real(KNDR),allocatable  ::  vertex_y(:,:)
      real(KNDR),allocatable  ::  area(:,:)
      real(KNDR),allocatable  ::  volm(:,:)
      
      character(lenofchar),allocatable  :: node_type_char(:,:)
      real(kndr),allocatable            :: node_type_direction(:,:,:)
      real(kndr),allocatable            :: node_type_velocity(:,:,:)
      real(kndr),allocatable            :: node_type_pressure(:,:)
         
   End Type Mesh_Cyl_2D
   
End Module COMMON_TYPE_Mesh


