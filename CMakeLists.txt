cmake_minimum_required(VERSION 3.5)

project(erfinv LANGUAGES C CXX)
set(CMAKE_CXX_STANDARD 11)

## Use Catch2 for unit testing
set(CATCH_VERSION "v2.12.2")
set(CATCH_SHA256 "a206bf2fef4059e5f805b448ddd12a29e41305c1ba02261286c2858bebdde9af")
file(DOWNLOAD
  https://raw.githubusercontent.com/catchorg/Catch2/${CATCH_VERSION}/single_include/catch2/catch.hpp
  ${CMAKE_BINARY_DIR}/external/catch/catch.hpp
  EXPECTED_HASH SHA256=${CATCH_SHA256}
)

add_library(catch INTERFACE)
target_include_directories(catch INTERFACE ${CMAKE_BINARY_DIR}/external)

## The erfinv library
add_library(erfinv erfinv.c)

## Unit tests
enable_testing()

add_executable(ppnd_test ppnd_test.cpp)
target_link_libraries(ppnd_test erfinv catch)
add_test(ppnd_test ppnd_test)

add_executable(erfinv_test erfinv_test.cpp)
target_link_libraries(erfinv_test erfinv catch)
add_test(erfinv_test erfinv_test)
