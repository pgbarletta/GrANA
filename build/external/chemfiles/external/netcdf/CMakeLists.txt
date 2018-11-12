cmake_minimum_required(VERSION 2.8.11)

project(netcdf-chemfiles C)
set(VERSION "4.4.1")

if(NOT WIN32)
    option(NETCDF_ENABLE_STDIO "Use stdio instead of POSIX IO (ex. on the Cray)" OFF)
    if(NETCDF_ENABLE_STDIO)
        set(USE_STDIO ON CACHE BOOL "")
    endif()

    option(NETCDF_ENABLE_FFIO "Use ffio instead of POSIX IO" OFF)
    if(NETCDF_ENABLE_FFIO)
        set(USE_FFIO ON CACHE BOOL "")
    endif()
endif()

# auto-configure style checks, other CMake modules.
include(CheckIncludeFile)
include(CheckTypeSize)
include(CheckFunctionExists)

# Library include checks
check_include_file("math.h"         HAVE_MATH_H)
check_include_file("unistd.h"       HAVE_UNISTD_H)
check_include_file("alloca.h"       HAVE_ALLOCA_H)
check_include_file("malloc.h"       HAVE_MALLOC_H)
check_include_file("ctype.h"        HAVE_CTYPE_H)
check_include_file("dirent.h"       HAVE_DIRENT_H)
check_include_file("dlfcn.h"        HAVE_DLFCN_H)
check_include_file("errno.h"        HAVE_ERRNO_H)
check_include_file("fcntl.h"        HAVE_FCNTL_H)
check_include_file("getopt.h"       HAVE_GETOPT_H)
check_include_file("stdbool.h"      HAVE_STDBOOL_H)
check_include_file("locale.h"       HAVE_LOCAL_H)
check_include_file("stdint.h"       HAVE_STDINT_H)
check_include_file("stdio.h"        HAVE_STDIO_H)
check_include_file("stdlib.h"       HAVE_STDLIB_H)
check_include_file("stdarg.h"       HAVE_STDARG_H)
check_include_file("strings.h"      HAVE_STRINGS_H)
check_include_file("signal.h"       HAVE_SIGNAL_H)
check_include_file("sys/dir.h"      HAVE_SYS_DIR_H)
check_include_file("sys/ndir.h"     HAVE_SYS_NDIR_H)
check_include_file("sys/param.h"    HAVE_SYS_PARAM_H)
check_include_file("sys/stat.h"     HAVE_SYS_STAT_H)
check_include_file("sys/time.h"     HAVE_SYS_TIME_H)
check_include_file("sys/types.h"    HAVE_SYS_TYPES_H)
check_include_file("sys/wait.h"     HAVE_SYS_WAIT_H)
check_include_file("sys/resource.h" HAVE_SYS_RESOURCE_H)
check_include_file("fcntl.h"        HAVE_FCNTL_H)
check_include_file("inttypes.h"     HAVE_INTTYPES_H)
check_include_file("pstdint.h"      HAVE_PSTDINT_H)
check_include_file("endian.h"       HAVE_ENDIAN_H)
check_include_file("BaseTsd.h"      HAVE_BASETSD_H)
check_include_file("stddef.h"       HAVE_STDDEF_H)

# Type checks
check_type_size("void*"     SIZEOF_VOIDSTAR)
check_type_size("char"      SIZEOF_CHAR)
check_type_size("double"    SIZEOF_DOUBLE)
check_type_size("float"     SIZEOF_FLOAT)
check_type_size("int"       SIZEOF_INT)
check_type_size("long"      SIZEOF_LONG)
check_type_size("long long" SIZEOF_LONG_LONG)
check_type_size("off_t"     SIZEOF_OFF_T)
check_type_size("off64_t"   SIZEOF_OFF64_T)
check_type_size("short"     SIZEOF_SHORT)
check_type_size("size_t"    SIZEOF_SIZE_T)
check_type_size("ssize_t"   SIZEOF_SSIZE_T)
# __int64 is used on Windows for large file support.
check_type_size("__int64"   SIZEOF___INT_64)
check_type_size("uchar"     SIZEOF_UCHAR)
check_type_size("int64_t"   SIZEOF_INT64_T)
check_type_size("uint64_t"  SIZEOF_UINT64_T)

# On windows systems, we redefine off_t as __int64
# to enable LFS. This is true on 32 and 64 bit system.s
# We must redefine SIZEOF_OFF_T to match.
IF(MSVC AND SIZEOF___INT_64)
  SET(SIZEOF_OFF_T  ${SIZEOF___INT_64})
ENDIF()

IF(SIZEOF_SSIZE_T)
  SET(HAVE_SSIZE_T 1)
ELSE()
  check_type_size("SSIZE_T" SIZEOF_SSIZE_T)
  IF(SIZEOF_SSIZE_T)
    SET(HAVE_SSIZE_T 1)
  ENDIF()
ENDIF()

# Check for various functions.
check_function_exists(fsync         HAVE_FSYNC)
check_function_exists(strlcat       HAVE_STRLCAT)
check_function_exists(strerror      HAVE_STRERROR)
check_function_exists(snprintf      HAVE_SNPRINTF)
check_function_exists(strchr        HAVE_STRCHR)
check_function_exists(strrchr       HAVE_STRRCHR)
check_function_exists(strcat        HAVE_STRCAT)
check_function_exists(strcpy        HAVE_STRCPY)
check_function_exists(strdup        HAVE_STRDUP)
check_function_exists(strndup       HAVE_STRNDUP)
check_function_exists(strcasecmp    HAVE_STRCASECMP)
check_function_exists(strtod        HAVE_STRTOD)
check_function_exists(strtoll       HAVE_STRTOLL)
check_function_exists(strtoull      HAVE_STROULL)
check_function_exists(strstr        HAVE_STRSTR)
check_function_exists(mkstemp       HAVE_MKSTEMP)
check_function_exists(rand          HAVE_RAND)
check_function_exists(random        HAVE_RANDOM)
check_function_exists(gettimeofday  HAVE_GETTIMEOFDAY)
check_function_exists(memmove       HAVE_MEMMOVE)
check_function_exists(getpagesize   HAVE_GETPAGESIZE)
check_function_exists(sysconf       HAVE_SYSCONF)
check_function_exists(getrlimit     HAVE_GETRLIMIT)
check_function_exists(_filelengthi64 HAVE_FILE_LENGTH_I64)

configure_file(
    config.h.in
    ${CMAKE_CURRENT_BINARY_DIR}/include/config.h
)

SET(NETCDF_SOURCES
    # netcdf entry point
    liblib/nc_initialize.c
    # netcdf3 sources
    libsrc/attr.c libsrc/dim.c libsrc/lookup3.c libsrc/nc_hashmap.c
    libsrc/nc3dispatch.c libsrc/nc3internal.c libsrc/ncio.c libsrc/ncx.c
    libsrc/putget.c libsrc/v1hpg.c libsrc/var.c
    # dispatch sources
    libdispatch/datt.c libdispatch/dattget.c libdispatch/dattinq.c
    libdispatch/dattput.c libdispatch/dcopy.c libdispatch/ddim.c
    libdispatch/ddispatch.c libdispatch/derror.c libdispatch/dfile.c
    libdispatch/dinternal.c libdispatch/dparallel.c libdispatch/dstring.c
    libdispatch/dutf8proc.c libdispatch/dvar.c libdispatch/dvarget.c
    libdispatch/dvarinq.c libdispatch/dvarput.c libdispatch/nc.c
    libdispatch/ncbytes.c libdispatch/nchashmap.c libdispatch/nclist.c
    libdispatch/nclistmgr.c libdispatch/nclog.c libdispatch/nctime.c
    libdispatch/ncuri.c
)

if(USE_FFIO)
    set(NETCDF_SOURCES ${NETCDF_SOURCES} libsrc/ffio.c)
elseif(USE_STDIO)
   set(NETCDF_SOURCES ${NETCDF_SOURCES} libsrc/ncstdio.c)
else()
    set(NETCDF_SOURCES ${NETCDF_SOURCES} libsrc/posixio.c)
endif()

add_library(netcdf OBJECT ${NETCDF_SOURCES})
target_include_directories(netcdf PUBLIC include ${CMAKE_CURRENT_BINARY_DIR}/include)
