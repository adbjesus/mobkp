# This file sets up GLPK for CMake. Once done this will define
#  GLPK_FOUND             - system has GLPK lib
#  GLPK_INCLUDE_DIR       - the GLPK include directory
#  GLPK_LIBRARIES         - Link these to use GLPK

if (NOT GLPK_FOUND)
  find_path(GLPK_INCLUDE_DIR NAMES glpk.h)
  find_library(GLPK_LIBRARIES NAMES libglpk glpk)
  if(GLPK_LIBRARIES AND GLPK_INCLUDE_DIR)
    set(GLPK_FOUND TRUE)
  endif()
endif()
