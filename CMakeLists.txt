cmake_minimum_required(VERSION 3.17)

if(NOT DEFINED PROJECT_NAME)
  set(MOBKP_NOT_SUBPROJECT ON)
endif()

project(MOBKP
  VERSION 0.1.0
  DESCRIPTION "Algorithms for the Multi-Objective Binary Knapsack Problem"
  LANGUAGES CXX
)


# Set warning flags
if (CMAKE_CXX_COMPILER_ID MATCHES "MSVC")
  set(MOBKP_CXX_WARN_FLAGS /W4)
else()
  # TODO check if these all work in clang, or if we want different ones
  set(MOBKP_CXX_WARN_FLAGS
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
    -Wvector-operation-performance
    -Wshadow
    -Wnull-dereference
    -Wduplicated-cond
  )
endif ()

# Project options
# By default build apps and documentation
option(MOBKP_BUILD_APPS "Build the solver." ON)
option(MOBKP_BUILD_DOCS "Build the documentation." ON)
# By default build tests and benchmarks only if not a subproject
option(MOBKP_BUILD_TESTS "Build the tests." ${MOBKP_NOT_SUBPROJECT})
option(MOBKP_BUILD_BENCHES "Build the benchmarks." ${MOBKP_NOT_SUBPROJECT})

# Extra stuff if this is not being used as a subproject
if(MOBKP_NOT_SUBPROJECT)
  # Turn off C++ extensions
  set(CMAKE_CXX_EXTENSIONS OFF)
  # Create compile commands json for clang tools
  set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
  # Include testing
  include(CTest)
endif()

# Define library target
add_library(mobkp INTERFACE)
target_compile_features(mobkp INTERFACE cxx_std_20)
target_include_directories(mobkp INTERFACE
  $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
  $<INSTALL_INTERFACE:include>
  )

find_package(mooutils QUIET)
if (NOT mooutils_FOUND)
  include(FetchContent)
  FetchContent_Declare(
    mooutils
    GIT_REPOSITORY git@git.adbjesus.com:mooutils.git
  )
  FetchContent_MakeAvailable(mooutils)
  add_library(mooutils::mooutils ALIAS mooutils)
endif()
target_link_libraries(mobkp INTERFACE mooutils::mooutils)

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
  # TODO add_subdirectory(benches)
endif()
