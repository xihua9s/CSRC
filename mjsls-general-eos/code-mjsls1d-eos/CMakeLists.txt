# 强制设置编译器 FC = ifort
#set(ENV{FC} "/opt/intel/composer_xe_2011_sp1.11.339/bin/ia32/ifort")
#set(ENV{FC} "/opt/intel/composer_xe_2011_sp1.11.339/bin/intel64/ifort")
#set(ENV{FC} "C:/Program Files (x86)/Intel/Composer XE 2013/bin/intel64/ifort")

CMAKE_MINIMUM_REQUIRED(VERSION 2.8...3.10.0) 

PROJECT(Xihua)
enable_language(Fortran)

# 设置编译类型 debug 或者 release
#SET(CMAKE_BUILE_TYPE DEBUG)
SET(CMAKE_BUILE_TYPE RELEASE)

# 设置编译器的类型
#SET(CMAKE_IFORT_FLAGS_DEBUG)

#include(CMAKEForceCompiler)
#CMAKE_FORCE_C_COMPILER(icc "INTEL C COMPILER")
#CMAKE_FORCE_CXX_COMPILER(ICPC "INTEL C++ COMPILER")
#CMAKE_FORCE_FORTRAN_COMPILER(IFORT "INTEL FORTRAN COMPILER")

set(CMAKE_Fortran_FLAGS " -O3 -ffast-math -Wuninitialized -fbounds-check")

#head file path
#INCLUDE_DIRECTORIES(
#   include
# )

#source directory
#AUX_SOURCE_DIRECTORY(src DIR_SRCS)

file(GLOB scr_files
	./src/main.f90
	./src/common/common_inicmt.f90
	./src/common/common_meshphys.f90
./src/posttreaty/output_mesh_phys.f90
./src/pretreaty/MDU_Initial.f90
./src/pretreaty/MDU_Initial_Kidder.f90
./src/pretreaty/MDU_Initial_lazarus.f90
./src/pretreaty/MDU_Initial_noh.f90
./src/pretreaty/MDU_Initial_s123.f90
./src/pretreaty/MDU_Initial_s321.f90
./src/pretreaty/MDU_Initial_sedov.f90
./src/pretreaty/MDU_Initial_JWL.f90
./src/pretreaty/MDU_Initial_sod.f90
./src/pretreaty/MDU_Initial_shuosher.f90
./src/pretreaty/MDU_Initial_twoblast.f90
./src/pretreaty/MDU_Initial_lax.f90
./src/pretreaty/MDU_Initial_leblanc.f90
./src/pretreaty/MDU_Initial_stiff.f90
./src/pretreaty/MDU_Initial_TNT.f90
./src/solver/MDU_exact_riemann_slover.f90
./src/solver/MDU_fluid_euler_eqn_solver.f90
./src/solver/cal_dt.f90
./src/solver/fluid_euler_eqn_basic.f90
./src/solver/centered/mdu_fluid_euler_eqn_solver_centered.f90
./src/solver/centered/mdu_fluid_euler_eqn_solver_centered_node_solver.f90
./src/solver/centered/mdu_fluid_euler_eqn_solver_centered_update.f90
./src/solver/highorder/highorder.f90
./src/solver/staggered/*.f90
./src/solver/eos/*.f90
)

#set environment variable
#SET(TEST_MATH
#${DIR_SRCS}
# )


#add executable file
ADD_EXECUTABLE(jos ${scr_files})

#add link library
#TARGET_LINK_LIBRARIES(../bin/bin m)
