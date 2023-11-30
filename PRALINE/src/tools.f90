

subroutine cal_length_normal(x1,x2,y1,y2,alength,anormalx,anormaly)
	implicit none
	real*8 :: x1,x2,y1,y2
	real*8 :: alength,anormalx,anormaly

	call Length(x1,y1,x2,y2,alength)
	call Normal(x1,y1,x2,y2,anormalx,anormaly)
	anormalx = -anormalx
	anormaly = -anormaly		

end subroutine cal_length_normal

subroutine Length(x1,y1,x2,y2,alength)
   implicit none
   real*8,intent(in) :: x1,y1,x2,y2
   real*8,intent(out) :: alength

   alength = dsqrt((x1-x2)**2+(y1-y2)**2)    
   
end subroutine Length


subroutine Normal(x1,y1,x2,y2,anormalx,anormaly)
   implicit none
   real*8,intent(in) :: x1,y1,x2,y2
   real*8,intent(out) :: anormalx,anormaly
   real*8 :: tempt

   tempt=sqrt((x1-x2)**2+(y1-y2)**2)
   anormalx=(y1-y2)/tempt
   anormaly=(x2-x1)/tempt
   
end subroutine Normal


subroutine cal_eos(eos,gamar)
   implicit none 
   integer,intent(in)  :: eos 
   real*8,intent(out)  :: gamar

   select case (eos)
   case (715)
      gamar = 7.d0 / 5.d0 
   case (513)
      gamar = 5.d0 / 3.d0 
   case (815)
      gamar = 8.d0 / 5.d0 
   case (312)
      gamar = 3.d0 / 2.d0
   case default
      stop 'eos,errorela'
   end select 
end subroutine cal_eos
      


