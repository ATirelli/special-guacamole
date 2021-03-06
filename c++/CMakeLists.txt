cmake_minimum_required(VERSION 3.15)
project(demux)

set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -Wall")
set(CMAKE_CXX_STANDARD 14)
find_package(Boost)
find_package(HDF5)

add_subdirectory(src)
#add_subdirectory(tst)
