cmake_minimum_required(VERSION 3.10.0)
project(GRANA VERSION 0.1.0 LANGUAGES CXX)

set(GRANA_CGAL_DIR "" CACHE STRING
"CGAL installation folder.")
set(GRANA_FMT_DIR "" CACHE STRING
"fmt installation folder.")
set(GRANA_CHEMFILES_DIR "" CACHE STRING
"Chemfiles installation folder. It should have an include and a lib folder.")

set(SOURCE_DIR "${CMAKE_INSTALL_PREFIX}/src")
set(INCLUDE_DIRECTORIES "${CMAKE_INSTALL_PREFIX}/include")

set(INCLUDES
    include/GRANA/utils.hpp
    include/GRANA/primitives.hpp
    include/GRANA/continuous.hpp
    include/GRANA/grid.hpp
include/GRANA/main.hpp)
set(SOURCES
    src/utils.cpp
    src/primitives.cpp
    src/continuous.cpp
    src/grid.cpp
src/main.cpp)

add_executable(GRANA ${SOURCES})

if(GRANA_CGAL_DIR STREQUAL "")
    message("Trying to automatically locate CGAL library.")
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

target_include_directories(GRANA
    PUBLIC
        $<INSTALL_INTERFACE:include>
    PRIVATE
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
        /home/pbarletta/chemfiles08/build/install/include
)
target_link_libraries(GRANA CGAL)
target_link_libraries(GRANA gmp)
target_link_libraries(GRANA mpfr)
target_link_libraries(GRANA fmt::fmt)
target_link_libraries(GRANA ${chemfiles})
target_link_libraries(GRANA netcdf)

target_compile_options(GRANA PRIVATE -Wall)
target_compile_features(GRANA PRIVATE cxx_std_17)
