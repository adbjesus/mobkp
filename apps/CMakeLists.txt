find_package(fmt REQUIRED)
find_package(Boost REQUIRED)

find_package(CLI11 QUIET)
if (NOT CLI11_FOUND)
  include(FetchContent)
  FetchContent_Declare(
    cli11
    GIT_REPOSITORY https://github.com/CLIUtils/CLI11
    GIT_TAG        v2.1.2
  )
  FetchContent_MakeAvailable(cli11)
endif()

add_executable(mobkp-bin mobkp.cpp)
target_link_libraries(mobkp-bin mobkp mooutils::mooutils fmt::fmt Boost::headers CLI11::CLI11)
target_compile_options(mobkp-bin PRIVATE ${MOBKP_CXX_WARN_FLAGS})
set_target_properties(mobkp-bin PROPERTIES OUTPUT_NAME mobkp)
install(TARGETS mobkp-bin)
