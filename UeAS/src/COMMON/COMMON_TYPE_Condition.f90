

Module COMMON_TYPE_Condition
   use COMMON_const_kind
   implicit none
   
   Type :: COMMON_TYPE_CDT
   
      character(36) :: Test_Case      
      character(36) :: Cat_Cyl        
      character(36) :: Gou_Xing    
      character(36) :: Origin_YesNo  
      character(36) :: Angle    
      integer(KNDI) :: nx       
      integer(KNDI) :: ny       
      integer(KNDI) :: It_num   
      real(KNDR)    :: Tstop   
      real(KNDR)    :: CFL     
      integer(KNDI) :: It_out     
      character(36) :: Node_Solver   
      character(36) :: CV_AW      
      integer(KNDI) :: EOS         
      integer(KNDI) :: it=0         
      real(KNDR)    :: time=0.d0   
      real(KNDR)    :: dt=0.d0    
      
      real(KNDR)    :: foucstime=0.d0      
   
   End Type COMMON_TYPE_CDT
   
End Module COMMON_TYPE_Condition
