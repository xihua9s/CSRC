
program main
   use MDU_Initial
   use MDU_fluid_euler_eqn_solver

   implicit none

   type(inicmt)   :: ini_cmt 
   type(meshphys) :: mesh_phys   

   print*, '---> Initial(ini_cmt,mesh_phys)'
   Call Initial(ini_cmt,mesh_phys)
   print*, '<--- Initial(ini_cmt,mesh_phys)'

   print*, '---> fluid_euler_eqn_solver(ini_cmt,mesh_phys)'
   Call fluid_euler_eqn_solver(ini_cmt,mesh_phys)

end program main
