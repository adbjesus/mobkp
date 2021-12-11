find_package(fmt REQUIRED)
find_package(CLI11 REQUIRED)
find_package(Boost REQUIRED)

add_executable(mobkp-bin mobkp.cpp)
target_link_libraries(mobkp-bin mobkp mooutils::mooutils fmt::fmt Boost::headers CLI11::CLI11)
target_compile_options(mobkp-bin PRIVATE ${MOBKP_CXX_WARN_FLAGS})
set_target_properties(mobkp-bin PROPERTIES OUTPUT_NAME mobkp)
install(TARGETS mobkp-bin)
