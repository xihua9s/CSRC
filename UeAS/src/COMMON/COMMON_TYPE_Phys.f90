

Module COMMON_TYPE_Phys
   use COMMON_const_kind
   implicit none
   
   Type :: Phys_Cyl_2D
   
      real(KNDR),allocatable  ::  den(:,:)
      real(KNDR),allocatable  ::  tau(:,:)
      real(KNDR),allocatable  ::  vel_x(:,:)
      real(KNDR),allocatable  ::  vel_y(:,:)
      real(KNDR),allocatable  ::  pre(:,:)
      real(KNDR),allocatable  ::  mass(:,:)
      real(KNDR),allocatable  ::  sound(:,:)
      real(KNDR),allocatable  ::  enin(:,:)
      real(KNDR),allocatable  ::  ener(:,:)
      integer(KNDI),allocatable  ::  EOS(:,:)
      
      real(KNDR),allocatable  ::  vstar_x(:,:)
      real(KNDR),allocatable  ::  vstar_y(:,:)
      real(KNDR),allocatable  ::  Pstar8(:,:,:)  
      real(KNDR),allocatable  ::  Pstar4(:,:,:)  
      real(KNDR),allocatable  ::  PstarC1(:,:)    
      real(KNDR),allocatable  ::  PstarC4(:,:,:) 
      real(KNDR),allocatable  ::  PstarC8(:,:,:) 
      
   End Type Phys_Cyl_2D
   
End Module COMMON_TYPE_Phys


