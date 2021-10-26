find_package(fmt REQUIRED)
find_package(CLI11 REQUIRED)
find_package(Boost REQUIRED)

add_executable(mobkpsolve mobkpsolve.cpp)
target_link_libraries(mobkpsolve mobkp mooutils::mooutils fmt::fmt Boost::headers)
target_compile_options(mobkpsolve PRIVATE ${MOBKP_CXX_WARN_FLAGS})
