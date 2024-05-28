# 强制设置编译器 FC = ifort
#set(ENV{FC} "/opt/intel/composer_xe_2011_sp1.11.339/bin/ia32/ifort")
#set(ENV{FC} "/opt/intel/composer_xe_2011_sp1.11.339/bin/intel64/ifort")
#set(ENV{FC} "C:/Program Files (x86)/Intel/Composer XE 2013/bin/intel64/ifort")

CMAKE_MINIMUM_REQUIRED(VERSION 2.8) 

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
	./src/*.f90
)

#set environment variable
#SET(TEST_MATH
#${DIR_SRCS}
# )


#add executable file
ADD_EXECUTABLE(jos ${scr_files})

#add link library
#TARGET_LINK_LIBRARIES(../bin/bin m)
