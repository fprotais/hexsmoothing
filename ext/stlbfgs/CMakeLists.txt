cmake_minimum_required(VERSION 3.21)
project(stlbfgs)
# modified CMakeLists

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)


set(RELATIVE_BIN_DIR bin/)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/${RELATIVE_BIN_DIR}/)
include_directories(${CMAKE_CURRENT_SOURCE_DIR})


file(GLOB SOURCES *.cpp *.h)
add_library(stlbfgs ${SOURCES})

set_target_properties(stlbfgs PROPERTIES SOVERSION 1.0)


