# This file sets up GLPK for CMake. Once done this will define
#  GLPK_FOUND             - system has GLPK lib
#  GLPK_INCLUDE_DIR       - the GLPK include directory
#  GLPK_LIBRARIES         - Link these to use GLPK

if (NOT mooutils_FOUND)
  include(FetchContent)
  FetchContent_Declare(
    mooutils
    GIT_REPOSITORY git@git.adbjesus.com:mooutils.git
  )
  FetchContent_MakeAvailable(mooutils)
  add_library(mooutils::mooutils ALIAS mooutils)
endif()
