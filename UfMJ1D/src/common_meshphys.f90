
module common_meshphys

	type meshphys
		real*8,allocatable  :: vertex(:)

		real*8,allocatable  :: con_phys(:,:)
		real*8,allocatable  :: pre(:)
		real*8,allocatable  :: sos(:)
		real*8,allocatable  :: ein(:)
		real*8,allocatable  :: vel_p(:)
		real*8,allocatable  :: mass(:)
		real*8,allocatable  :: pstar(:)
		real*8,allocatable  :: vol_size_t0(:)
		real*8,allocatable  :: vol_size_old(:)
		real*8,allocatable  :: vol_size_new(:)
		real*8,allocatable  :: entropy(:)
		real*8,allocatable  :: volnmin1(:)
		real*8,allocatable  :: voln_ratio(:)
		real*8,allocatable  :: den0(:)
		real*8,allocatable  :: vol0(:)
		real*8,allocatable  :: vol(:)
		real*8,allocatable  :: tau(:)
		real*8,allocatable  :: tau_old(:)
		real*8,allocatable  :: q(:)

		integer,allocatable :: eos(:)
		integer :: bnd_con(2)     
		real*8  :: bnd_vel_pre(4) 

	endtype meshphys

end module common_meshphys

