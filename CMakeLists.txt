cmake_minimum_required(VERSION 3.1)
project(hexsmoothing)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

include_directories(ext/stlbfgs)

add_subdirectory(ext/stlbfgs)
add_subdirectory(ext/ultimaille)



find_package(OpenMP)
if(OPENMP_FOUND)
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
  set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

if (WIN32) 
	set(BUILD_SHARED_LIBS OFF)
endif()

file(GLOB SOURCES utils/*.cpp utils/*.h )
add_library(libBasicUtilies ${SOURCES})

file(GLOB SOURCES lib/*.cpp lib/*.h )
add_library(libhexsmoothing ${SOURCES})
target_link_libraries(libhexsmoothing libBasicUtilies stlbfgs)
target_include_directories(libhexsmoothing PUBLIC ${CMAKE_CURRENT_SOURCE_DIR} ext/)
target_include_directories(libhexsmoothing PUBLIC ext/)

add_executable(innerSmoother bin/innerSmoothing.cpp)
target_link_libraries(innerSmoother libhexsmoothing ultimaille)
target_include_directories(innerSmoother PUBLIC ext/ultimaille)

add_executable(mixedSmoothing bin/mixedSmoothing.cpp)
target_link_libraries(mixedSmoothing libhexsmoothing ultimaille)
target_include_directories(mixedSmoothing PUBLIC ext/ultimaille)

if(OpenMP_CXX_FOUND)
    target_link_libraries(libhexsmoothing OpenMP::OpenMP_CXX)
endif()

