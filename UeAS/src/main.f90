
program main

   use COMMON_TYPE_Condition
   implicit none
   
   Type(COMMON_TYPE_CDT) :: Ini_CDT  
   
   Call Initial_Condition(Ini_CDT)

   Call Test_Benchmarck_Case(Ini_CDT)

end program main



