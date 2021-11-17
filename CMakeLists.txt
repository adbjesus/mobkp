cmake_minimum_required(VERSION 3.21)

if(NOT DEFINED PROJECT_NAME)
  set(MOBKP_MASTER_PROJECT ON)
else()
  set(MOBKP_MASTER_PROJECT OFF)
endif()

project(
  MOBKP
  VERSION 0.1.0
  DESCRIPTION "Anytime algorithms for the Multi-Objective Binary Knapsack Problem"
  LANGUAGES CXX
)

# Project options
option(MOBKP_DEVELOPMENT_BUILD "Development build" OFF)
option(MOBKP_BUILD_APPS "Build the apps" ${MOBKP_MASTER_PROJECT})
option(MOBKP_BUILD_DOCS "Build the documentation" ${MOBKP_MASTER_PROJECT})
option(MOBKP_BUILD_TESTS "Build tests" ${MOBKP_DEVELOPMENT_BUILD})
option(MOBKP_BUILD_BENCHMARKS "Build benchmarks" ${MOBKP_DEVELOPMENT_BUILD})
option(MOBKP_ENABLE_WARNINGS "Enable warnings" ${MOBKP_DEVELOPMENT_BUILD})
option(MOBKP_ENABLE_WERROR "Enable -Werror or /WX" ${MOBKP_DEVELOPMENT_BUILD})

# Warning flags
set(MOBKP_CXX_WARN_FLAGS "")
if (MOBKP_ENABLE_WARNINGS)
  if (CMAKE_CXX_COMPILER_ID MATCHES "MSVC")
    list(APPEND MOBKP_CXX_WARN_FLAGS /W4)
  elseif(CMAKE_CXX_COMPILER_ID MATCHES "GNU" OR CMAKE_CXX_COMPILER_ID MATCHES "Clang")
    list(APPEND MOBKP_CXX_WARN_FLAGS
      -Wall
      -Wextra
      -Wpedantic
      -pedantic-errors
      -Wold-style-cast
      -Wundef
      -Wredundant-decls
      -Wcast-qual
      -Wmissing-include-dirs
      -Wdisabled-optimization
      -Wconversion
      -Wdouble-promotion
      -Wshadow
      -Wnull-dereference
    )
    if(CMAKE_CXX_COMPILER_ID MATCHES "GNU")
      list(APPEND MOBKP_CXX_WARN_FLAGS
        -Wvector-operation-performance
        -Wduplicated-cond
      )
    endif()
  else()
    message(WARNING "Unknown compiler (not GNU, Clang, or MSVC), not setting any warning flags.")
    set(MOBKP_CXX_WARN_FLAGS "")
  endif()
endif()

if(MOBKP_ENABLE_WERROR)
  list(APPEND MOBKP_CXX_WARN_FLAGS
    $<$<CXX_COMPILER_ID:GNU,Clang>:-Werror>
    $<$<CXX_COMPILER_ID:MSVC>:/WX>
  )
endif()

# Extra stuff if this is the master project
if(MOBKP_MASTER_PROJECT)
  # Turn off C++ extensions
  set(CMAKE_CXX_EXTENSIONS OFF)
  # Create compile commands json for clang tools
  set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
  # Include CTest
  include(CTest)
endif()

# Find dependencies
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}/cmake/modules")

find_package(mooutils QUIET)
# If mooutils not found fallback to github
if (NOT mooutils_FOUND)
  include(FetchContent)
  FetchContent_Declare(
    mooutils
    GIT_REPOSITORY git@git.adbjesus.com:mooutils.git
  )
  FetchContent_MakeAvailable(mooutils)
  add_library(mooutils::mooutils ALIAS mooutils)
endif()

find_package(glpk REQUIRED)
find_package(Boost REQUIRED)

include(GNUInstallDirs)

# Define library target
add_library(mobkp INTERFACE)
target_compile_features(mobkp INTERFACE cxx_std_20)
target_include_directories(mobkp INTERFACE
  $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
  $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
  )
target_link_libraries(mobkp INTERFACE mooutils::mooutils glpk Boost::headers)

# The library is an install target
install(TARGETS mobkp)
install(DIRECTORY "${PROJECT_SOURCE_DIR}/include/mobkp" DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

# PkgConfig file
configure_file(
  "${PROJECT_SOURCE_DIR}/cmake/mobkp.pc.in"
  "${PROJECT_BINARY_DIR}/mobkp.pc"
  @ONLY
)
install(
  FILES "${PROJECT_BINARY_DIR}/mobkp.pc"
  DESTINATION "${CMAKE_INSTALL_LIBDIR}/pkgconfig"
)

# Export
install(
  TARGETS mobkp
  EXPORT mobkp-targets
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
  INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
)

install(
  EXPORT mobkp-targets
  NAMESPACE mobkp::
  FILE mobkp-targets.cmake
  DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/mobkp"
)

include(CMakePackageConfigHelpers)
configure_package_config_file(
  "${PROJECT_SOURCE_DIR}/cmake/mobkp-config.cmake.in"
  "${PROJECT_BINARY_DIR}/mobkp-config.cmake"
  INSTALL_DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/mobkp"
)
write_basic_package_version_file(
  "${PROJECT_BINARY_DIR}/mobkp-config-version.cmake"
  VERSION "${PROJECT_VERSION}"
  COMPATIBILITY SameMinorVersion
)

install(
  FILES "${PROJECT_BINARY_DIR}/mobkp-config.cmake"
        "${PROJECT_BINARY_DIR}/mobkp-config-version.cmake"
  DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/mobkp"
)

# Define solver target
if(MOBKP_BUILD_APPS)
  add_subdirectory(apps)
endif()

# Define docs target
if(MOBKP_BUILD_DOCS)
  # add_subdirectory(docs)
endif()

# Define tests target
if(MOBKP_BUILD_TESTS)
  # add_subdirectory(tests)
endif()

# Define tests target
if(MOBKP_BUILD_BENCHES)
  # add_subdirectory(benches)
endif()
