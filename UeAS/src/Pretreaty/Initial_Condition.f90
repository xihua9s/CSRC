

Subroutine Initial_Condition(Ini_CDT)
   use COMMON_TYPE_Condition
   use COMMON_XXX_Constant
   implicit none
   
   Type(COMMON_TYPE_CDT),intent(out) :: Ini_CDT
   integer(KNDI) :: date_time(8)
   character(10) :: b(3)
   
   print*, '---> Initial_Condition(Ini_CDT)'
   
   open(1,file=FoldInput//'/Ini_CDT.ini')
   open(2,file=FoldOutput//'/Ini_CDT.ini')

   read(1,*)  Ini_CDT%Test_Case     
   print*,    Ini_CDT%Test_Case,'/Test Case/'  
   write(2,*) Ini_CDT%Test_Case,'/Test Case/'  

   read(1,*)  Ini_CDT%Cat_Cyl      
   print*,    Ini_CDT%Cat_Cyl,'/Cat_Cyl/'  
   write(2,*) Ini_CDT%Cat_Cyl,'/Cat_Cyl/'  

   read(1,*)  Ini_CDT%Gou_Xing      
   print*,    Ini_CDT%Gou_Xing,'/Gou_Xing/'  
   write(2,*) Ini_CDT%Gou_Xing,'/Gou_Xing/'  

   read(1,*)  Ini_CDT%Origin_YesNo      
   print*,    Ini_CDT%Origin_YesNo,'/Origin_YesNo/'  
   write(2,*) Ini_CDT%Origin_YesNo,'/Origin_YesNo/'  

   read(1,*)  Ini_CDT%Angle     
   print*,    Ini_CDT%Angle,'/Angle/'  
   write(2,*) Ini_CDT%Angle,'/Angle/'  

   read(1,*)  Ini_CDT%nx       
   print*,    Ini_CDT%nx, '/nx/'
   write(2,*) Ini_CDT%nx, '/nx/'

   read(1,*)  Ini_CDT%ny     
   print*,    Ini_CDT%ny, '/ny/'
   write(2,*) Ini_CDT%ny, '/ny/'

   read(1,*)  Ini_CDT%Tstop       
   print*,    Ini_CDT%Tstop, '/Tstop/'
   write(2,*) Ini_CDT%Tstop, '/Tstop/'


   read(1,*)  Ini_CDT%CFL       
   print*,   Ini_CDT%CFL,'/CFL/'
   write(2,*) Ini_CDT%CFL,'/CFL/'

   read(1,*)  Ini_CDT%It_num      
   print*,    Ini_CDT%It_num,'/Iterative Num/'
   write(2,*) Ini_CDT%It_num,'/Iterative Num/'

   read(1,*)  Ini_CDT%It_out              
   print*,   Ini_CDT%It_out,' /Iterative Out/' 
   write(2,*) Ini_CDT%It_out,' /Iterative Out/'
   
   read(1,*)  Ini_CDT%Node_Solver                
   print*,    Ini_CDT%Node_Solver,' /Node_Solver/' 
   write(2,*) Ini_CDT%Node_Solver,' /Node_Solver/'
   
   read(1,*)  Ini_CDT%CV_AW           
   print*,    Ini_CDT%CV_AW,' /CV_AW/'  
   write(2,*) Ini_CDT%CV_AW,' /CV_AW/' 
   
   read(1,*)  Ini_CDT%EOS          
   print*,    Ini_CDT%EOS,' /EOS/'  
   write(2,*) Ini_CDT%EOS,' /EOS/'    
   
   close(1)
   
   Call date_and_time(b(1),b(2),b(3),date_time)
   write(2,*)
   write(2,*) '==============================='
   write(2,*) 'Date:',b(1),b(2)
   write(2,*) '@copyright Xihua'
   write(2,*) '==============================='
   close(2)

   print*, '<--- Initial_Condition(Ini_CDT)'   
End Subroutine Initial_Condition