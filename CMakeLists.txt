cmake_minimum_required(VERSION 3.10.0)
project(grana VERSION 0.1.0 LANGUAGES CXX)

###################################################################
if(NOT WIN32)
  string(ASCII 27 Esc)
  set(ColourReset "${Esc}[m")
  set(ColourBold  "${Esc}[1m")
  set(Red         "${Esc}[31m")
  set(Green       "${Esc}[32m")
  set(Yellow      "${Esc}[33m")
  set(Blue        "${Esc}[34m")
  set(Magenta     "${Esc}[35m")
  set(Cyan        "${Esc}[36m")
  set(White       "${Esc}[37m")
  set(BoldRed     "${Esc}[1;31m")
  set(BoldGreen   "${Esc}[1;32m")
  set(BoldYellow  "${Esc}[1;33m")
  set(BoldBlue    "${Esc}[1;34m")
  set(BoldMagenta "${Esc}[1;35m")
  set(BoldCyan    "${Esc}[1;36m")
  set(BoldWhite   "${Esc}[1;37m")
endif()
###################################################################

set(GRANA_CGAL_DIR "" CACHE STRING
"CGAL installation folder.")
set(GRANA_FMT_DIR "" CACHE STRING
"fmt installation folder.")
set(GRANA_CHEMFILES_DIR "" CACHE STRING
"Chemfiles installation folder. It should have an include and a lib folder.")

set(SOURCE_DIR "${CMAKE_INSTALL_PREFIX}/src")
set(INCLUDE_DIRECTORIES "${CMAKE_INSTALL_PREFIX}/include")

set(INCLUDES
    include/grana/utils.hpp
    include/grana/primitives.hpp
    include/grana/continuous.hpp
    include/grana/grid.hpp
include/grana/main.hpp)
set(SOURCES
    src/utils.cpp
    src/primitives.cpp
    src/continuous.cpp
    src/grid.cpp
src/main.cpp)

add_executable(grana ${SOURCES})

if(GRANA_CGAL_DIR STREQUAL "")
    message("${Magenta}Trying to automatically locate CGAL library.")
    find_package(CGAL 4.5 CONFIG REQUIRED)
else()
    message("Not yet supported.")
endif()

if(GRANA_FMT_DIR STREQUAL "")
    message("Trying to automatically locate fmt library.")
    find_package(fmt 5.0 CONFIG REQUIRED)
else()
    message("Not yet supported.")
endif()

if(GRANA_CHEMFILES_DIR STREQUAL "")
    message("Trying to automatically locate Chemfiles library.")
    find_package(chemfiles 0.8 CONFIG REQUIRED)
else()
    message("Looking for a Chemfiles installation in ${GRANA_CHEMFILES_DIR}")
    set(CMAKE_PREFIX_PATH ${GRANA_CHEMFILES_DIR})
    find_library(chemfiles
        NAMES chemfiles
        HINT "${GRANA_CHEMFILES_DIR}/lib"
    )
    if(NOT chemfiles)
        message(FATAL_ERROR "Chemfiles was not found in the provided path.")
    else()
        message("Found Chemfiles!")
    endif()
endif()

message(${ColourReset})

target_include_directories(grana
    PUBLIC
        $<INSTALL_INTERFACE:include>
    PRIVATE
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
        /home/pbarletta/chemfiles08/build/install/include
)
target_link_libraries(grana CGAL)
target_link_libraries(grana gmp)
target_link_libraries(grana mpfr)
target_link_libraries(grana fmt::fmt)
target_link_libraries(grana ${chemfiles})
target_link_libraries(grana netcdf)

target_compile_options(grana PRIVATE -Wall)
target_compile_features(grana PRIVATE cxx_std_17)
