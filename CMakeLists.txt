cmake_minimum_required(VERSION 3.10.0)
project(grana VERSION 0.1.0 LANGUAGES CXX)

###################################################################
if(NOT WIN32)
  string(ASCII 27 Esc)
  set(BoldRed     "${Esc}[1;31m")
  set(BoldGreen   "${Esc}[1;32m")
  set(BoldYellow  "${Esc}[1;33m")
  set(BoldBlue    "${Esc}[1;34m")
  set(BoldMagenta "${Esc}[1;35m")
  set(BoldCyan    "${Esc}[1;36m")
  set(BoldWhite   "${Esc}[1;37m")
endif()
###################################################################

set(GRANA_CGAL_DIR "" CACHE STRING "CGAL installation folder.")

add_executable(grana "")
target_sources(grana
    PRIVATE
        src/main.cpp
        src/utils.cpp
        src/continuous_base.cpp
        src/continuous.cpp
        src/grid.cpp
        src/grid_base.cpp
        src/octree.cpp
        src/PDB.cpp
        include/GrANA/utils.hpp
        include/GrANA/continuous_primitives.hpp
        include/GrANA/continuous_base.hpp
        include/GrANA/continuous.hpp
        include/GrANA/grid_primitives.hpp
        include/GrANA/grid_base.hpp
        include/GrANA/grid.hpp
        include/GrANA/octree.hpp
        include/GrANA/PDB.hpp
)

message("${Magenta}Trying to automatically locate CGAL, gmp and mpfr libraries.${ColourReset}")
find_package(CGAL CONFIG REQUIRED)
target_link_libraries(grana PRIVATE CGAL)
target_link_libraries(grana PRIVATE gmp)
target_link_libraries(grana PRIVATE mpfr)

message("${Magenta}Building own chemfiles library.${ColourReset}")
add_subdirectory(external/chemfiles EXCLUDE_FROM_ALL)
target_link_libraries(grana PRIVATE chemfiles)

message("${Magenta}Using chemfiles's fmt header-only library.${ColourReset}")
include_directories(build/external/chemfiles/external/fmt/include)

# find_package(fmt)
# target_link_libraries(<your-target> fmt::fmt)

target_include_directories(grana
    PRIVATE
        $<INSTALL_INTERFACE:include>
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
)

target_compile_options(grana PRIVATE -Wall -Wextra)
target_compile_features(grana PRIVATE cxx_std_17)
