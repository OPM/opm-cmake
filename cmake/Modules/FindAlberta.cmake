set(ALBERTA_MAX_WORLD_DIM "3" CACHE STRING "Maximal world dimension to check for Alberta library.")
find_library(ALBERTA_LTDL_LIB
  NAMES ltdl
  PATH_SUFFIXES lib lib32 lib64
)
find_path(ALBERTA_INCLUDE_DIR
  NAMES alberta/alberta.h
  PATHS ${ALBERTA_ROOT}
  PATH_SUFFIXES alberta include NO_DEFAULT_PATH
  DOC "Include path of Alberta")
find_path(ALBERTA_INCLUDE_DIR
  NAMES
  alberta/alberta.h
  PATHS /usr/local /opt
  PATH_SUFFIXES alberta)
#look for libraries
find_library(ALBERTA_UTIL_LIB
  NAMES alberta_util alberta_utilities
  PATHS ${ALBERTA_ROOT}
  PATH_SUFFIXES lib lib32 lib64
  NO_DEFAULT_PATH)
find_library(ALBERTA_UTIL_LIB
  NAMES alberta_util alberta_utilities
  PATH_SUFFIXES lib lib32 lib64)

foreach(dim RANGE 1 ${ALBERTA_MAX_WORLD_DIM})
  find_library(ALBERTA_${dim}D_LIB alberta_${dim}d
    PATHS ${ALBERTA_ROOT}
    PATH_SUFFIXES lib lib32 lib64
    Cache FILEPATH DOC "Alberta lib for ${dim}D" NO_DEFAULT_PATH)
  find_library(ALBERTA_${dim}D_LIB alberta_${dim}d  PATH_SUFFIXES lib lib32 lib64)
  if(ALBERTA_${dim}D_LIB)
    set(ALBERTA_LIBRARIES ${ALBERTA_LIBRARIES} ${ALBERTA_${dim}D_LIB})
  endif()
endforeach(dim RANGE 1 9)

if(ALBERTA_LIBRARIES AND ALBERTA_INCLUDE_DIR)
  set(ALBERTA_INCLUDE_DIRS ${ALBERTA_INCLUDE_DIR})
  set(ALBERTA_LIBRARIES ${ALBERTA_LIBRARIES} ${ALBERTA_UTIL_LIB} ${ALBERTA_LTDL_LIB})
  set(ALBERTA_FOUND ON)
  set(Alberta_FOUND ON)
  set(HAVE_ALBERTA 1)
  set(DUNE_ALBERTA_VERSION 0x300)
else()
  set(ALBERTA_FOUND OFF)
  set(Alberta_FOUND OFF)
endif()

# PkgConfig targets are required if OPM modules are used as DUNE modules
# (e.g. in Debian's Packaging automatic tests)
find_package(PkgConfig)
if(PkgConfig_FOUND)
  set(_opm_alberta_bkup_path  ${CMAKE_PREFIX_PATH})
  set(CMAKE_PREFIX_PATH ${CMAKE_PREFIX_PATH} ${ALBERTA_ROOT})

  foreach(dim RANGE 1 ${ALBERTA_MAX_WORLD_DIM})
    pkg_check_modules(Alberta${dim}d IMPORTED_TARGET GLOBAL "alberta-grid_${dim}d")
  endforeach()

  set(CMAKE_PREFIX_PATH ${_opm_alberta_bkup_path})
endif()
