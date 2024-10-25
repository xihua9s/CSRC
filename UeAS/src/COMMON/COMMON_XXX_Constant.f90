
   
module COMMON_XXX_Constant
   use COMMON_const_kind
   implicit none
   
   real(KNDR),parameter :: eps = EPSILON(0.0_KNDR)
   real(KNDR),parameter :: lageps = -Huge(0.0_KNDR)
   real(KNDR),parameter :: ILageps = -999999999
   real(KNDR),parameter :: PI = 3.14159265359_KNDR
   
   character(9),parameter  :: FoldInput='InputData'
   character(10),parameter :: FoldOutput='OutputData'
   
end module COMMON_XXX_Constant
   
   


