CMAKE_MINIMUM_REQUIRED(VERSION 2.8) 

PROJECT(test_math)
enable_language(Fortran)


#head file path
#INCLUDE_DIRECTORIES(
#   include
# )

#source directory
#AUX_SOURCE_DIRECTORY(src DIR_SRCS)

file(GLOB scr_files
./src/COMMON/COMMON_TYPE_Condition.f90
./src/COMMON/COMMON_TYPE_EOS.f90
./src/COMMON/COMMON_TYPE_Mesh.f90
./src/COMMON/COMMON_TYPE_Phys.f90
./src/COMMON/COMMON_XXX_Constant.f90
./src/COMMON/COMMON_const_kind.f90
./src/Posttreaty/output.f90
./src/Pretreaty/Initial_Condition.f90
./src/Pretreaty/Initial_Mesh_Phys.f90
./src/Pretreaty/mdu_sod_2D_cat_cyl.f90
./src/Solver/Node_Solver.f90
./src/Solver/Solve_EulerEqn.f90
./src/Solver/Update_Mesh_Phys.f90
./src/Solver/UeAS/Node_Solver_UeAS.f90
./src/Solver/UeAS/Node_Solver_UeAS_BndyNode.f90
#./src/Solver/UeAS/Node_Solver_UeAS_BndyNode_Aspherical.f90
#./src/Solver/UeAS/Node_Solver_UeAS_BndyNode_Coggeshall.f90
#./src/Solver/UeAS/Node_Solver_UeAS_BndyNode_Kidder.f90
#./src/Solver/UeAS/Node_Solver_UeAS_BndyNode_Lazarus.f90
#./src/Solver/UeAS/Node_Solver_UeAS_BndyNode_Noh.f90
#./src/Solver/UeAS/Node_Solver_UeAS_BndyNode_Saltzmann.f90
#./src/Solver/UeAS/Node_Solver_UeAS_BndyNode_Sedov.f90
#./src/Solver/UeAS/Node_Solver_UeAS_BndyNode_Sod.f90
./src/Solver/UeAS/Update_Mesh_Phys_UeAS.f90
./src/Solver/get_dt.f90
./src/example/Test_Benchmarck_Case.f90
./src/main.f90
./src/tools/tools.f90
)

#set environment variable
#SET(TEST_MATH
#${DIR_SRCS}
# )


#add executable file
ADD_EXECUTABLE(jos ${scr_files})

#add link library
#TARGET_LINK_LIBRARIES(../bin/bin m)
# 
