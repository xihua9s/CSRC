
CMAKE_MINIMUM_REQUIRED(VERSION 2.8) 

PROJECT(Xihua)
enable_language(Fortran)

SET(CMAKE_BUILE_TYPE RELEASE)

set(CMAKE_Fortran_FLAGS " -O3 -ffast-math -Wuninitialized -fbounds-check")

file(GLOB scr_files
	./src/*.f90
)

ADD_EXECUTABLE(PRALINE ${scr_files})

