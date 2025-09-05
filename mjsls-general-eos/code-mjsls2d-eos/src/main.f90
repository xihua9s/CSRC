
!> 20240523-531 用面格式复现文献结果

!> 20250702 Plane 
!> 理解代码，找到粘性的计算方式
Program main
	use common_inicmt
	use mdu_mjsls_ca
	use mdu_mjsls_cy

	type(inicmt) :: ini_cmt

	call initial_cmt(ini_cmt)
	
	call mjsls_ca(ini_cmt,ini_cmt%nx,ini_cmt%ny)

	! call mjsls_cy_aw()
   
end program



