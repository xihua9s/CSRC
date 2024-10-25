

Module COMMON_TYPE_EOS
   use COMMON_const_kind
   implicit none
   
   real(KNDR)  ::  gammar
   
contains

Subroutine General_EOS(num)
   implicit none
   integer(KNDI),intent(in)  :: num
   
   Select Case(num)
   Case (1)
      gammar = 7.d0/5.d0
   Case (2)
      gammar = 5.d0/3.d0
   Case (3)
      gammar = 2.d0   
   Case default
      print*, 'eos num is:' , num
      stop 'The input num is ERROR, pls Call Xihua.'
   
   End Select

End subroutine General_EOS

End Module COMMON_TYPE_EOS


