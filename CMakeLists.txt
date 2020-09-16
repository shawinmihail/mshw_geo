cmake_minimum_required(VERSION 3.1)

project(mshw_geo)

include_directories(include)

add_library(mshw_geolib src/mshw_geolib.cpp)

add_executable(mshw_geo src/mshw_geolib_test.cpp)

target_link_libraries(mshw_geo mshw_geolib)
