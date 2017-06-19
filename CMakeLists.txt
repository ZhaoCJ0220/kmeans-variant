cmake_minimum_required(VERSION 3.5)

project(intellid_recolor)

SET(PACKAGE_VERSION "0.1.0")
SET(RECOLOR_MAJOR_VERSION 0)  # 0 -
SET(RECOLOR_MINOR_VERSION 1)  # 0 - 9
SET(RECOLOR_PATCH_VERSION 0)  # 0 - 99

MATH(EXPR RECOLOR_VERSION "${RECOLOR_MAJOR_VERSION} * 1000
     + ${RECOLOR_MINOR_VERSION} * 100 + ${RECOLOR_PATCH_VERSION}")

set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -std=c99 -pedantic -Wall -pthread")
set(CMAKE_C_FLAGS_DEBUG "-g -O0")
set(CMAKE_C_FLAGS_RELEASE "-O3")

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -pedantic -Wall -pthread")
set(CMAKE_CXX_FLAGS_DEBUG "-g -O0")
set(CMAKE_CXX_FLAGS_RELEASE "-O3")
#set(CMAKE_VERBOSE_MAKEFILE 1)

## Install the executable to <install_directory>/bin
set(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR}/bin)
## Install libraries (if any) to <install_directory>/bin
set(LIBRARY_OUTPUT_PATH ${PROJECT_BINARY_DIR}/pkg)

add_subdirectory(${PROJECT_SOURCE_DIR}/src)
