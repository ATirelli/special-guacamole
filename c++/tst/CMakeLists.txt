set(BINARY ${CMAKE_PROJECT_NAME}_tst)

file(GLOB_RECURSE TEST_SOURCES LIST_DIRECTORIES false *.h *.cpp)

set(SOURCES ${TEST_SOURCES})
find_package(Boost)
find_package(HDF5)
include_directories(${HDF5_INCLUDE_DIRS} ${Boost_INCLUDE_DIRS})

add_executable(${BINARY} ${TEST_SOURCES} test_signal_processing.cpp test_dtw.cpp)

add_test(NAME ${BINARY} COMMAND ${BINARY})

target_link_libraries(${BINARY} PUBLIC ${CMAKE_PROJECT_NAME}_lib gtest gmock
                      ${HDF5_LIBRARIES} ${Boost_LIBRARIES})