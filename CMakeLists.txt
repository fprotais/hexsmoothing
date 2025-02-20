cmake_minimum_required(VERSION 3.27)
project(hexsmoothing)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

# set(CMAKE_CXX_FLAGS "-Wall -Wextra")
set(CMAKE_CXX_FLAGS_DEBUG "-g")
set(CMAKE_CXX_FLAGS_RELEASE "-O3")

include_directories(ext/stlbfgs)
include_directories(ext/LBFGS-Lite/include)

add_subdirectory(ext/stlbfgs)
add_subdirectory(ext/ultimaille)

file(GLOB SOURCES utils/*.cpp utils/*.h )
add_library(libBasicUtilies ${SOURCES})

file(GLOB SOURCES lib/*.cpp lib/*.h )
add_library(libhexsmoothing ${SOURCES})
target_link_libraries(libhexsmoothing ultimaille libBasicUtilies stlbfgs)
target_include_directories(libhexsmoothing PUBLIC ${CMAKE_CURRENT_SOURCE_DIR} ext/)
target_include_directories(libhexsmoothing PUBLIC ext/ultimaille)

add_executable(innerSmoother bin/innerSmoothing.cpp)
target_link_libraries(innerSmoother libhexsmoothing ultimaille)
target_include_directories(innerSmoother PUBLIC ext/ultimaille)

add_executable(mixedSmoothing bin/mixedSmoothing.cpp)
target_link_libraries(mixedSmoothing libhexsmoothing ultimaille)
target_include_directories(mixedSmoothing PUBLIC ext/ultimaille)

option(PARALLEL_MIXED_SMOOTHING "" ON)

if (${PARALLEL_MIXED_SMOOTHING})
  find_package(OpenMP 4.0 REQUIRED)
  find_package(Eigen3 3.1.0 REQUIRED)

  file(GLOB SOURCES parallel_lib/*.cpp parallel_lib/*.h )
  add_library(libparallel ${SOURCES})
  target_link_libraries(libparallel libBasicUtilies libhexsmoothing OpenMP::OpenMP_CXX Eigen3::Eigen)

  add_executable(parallelMixedSmoothing bin/parallelMixedSmoothing.cpp)
  target_link_libraries(parallelMixedSmoothing libparallel)
endif()
