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
    message("${Magenta}Trying to automatically locate CGAL library.${ColourReset}")
    find_package(CGAL CONFIG REQUIRED)
    target_link_libraries(grana PRIVATE CGAL)
    target_link_libraries(grana PRIVATE gmp)
    target_link_libraries(grana PRIVATE mpfr)

    if(CGAL_FOUND)
        message("${Magenta}Found CGAL.${ColourReset}")
    elseif(NOT CGAL_FOUND)
        message("${Magenta}Could not found CGAL.${ColourReset}")
    endif()
else()
    message("${Magenta}Not yet supported.${ColourReset}")
endif()


if(GRANA_CHEMFILES_DIR STREQUAL "")
    message("${Magenta}Building own chemfiles library.${ColourReset}")
    add_subdirectory(external/chemfiles EXCLUDE_FROM_ALL)
    target_link_libraries(grana PRIVATE chemfiles)
else()
    message("${Magenta}Looking for a chemfiles installation in ${GRANA_CHEMFILES_DIR}.${ColourReset}")
    set(CMAKE_PREFIX_PATH ${GRANA_CHEMFILES_DIR})
    find_library(chemfiles
        NAMES chemfiles
        HINT "${GRANA_CHEMFILES_DIR}/lib"
    )
    if(NOT chemfiles)
        message(FATAL_ERROR "${Magenta}Chemfiles was not found in the provided path.${ColourReset}")
    else()
        message("${Magenta}Found Chemfiles!${ColourReset}")
    endif()
    target_link_libraries(grana PRIVATE ${chemfiles})
endif()

if(GRANA_FMT_DIR STREQUAL "")
    message("${Magenta}Using chemfiles's fmt header-only library.${ColourReset}")
    include_directories(build/external/chemfiles/external/fmt)
else()
    message("${Magenta}Not yet supported.${ColourReset}")
endif()

target_include_directories(grana
    PUBLIC
        $<INSTALL_INTERFACE:include>
    PRIVATE
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
        /home/pbarletta/chemfiles08/build/install/include
)

target_compile_options(grana PRIVATE -Wall)
target_compile_features(grana PRIVATE cxx_std_17)
