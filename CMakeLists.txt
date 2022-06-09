cmake_minimum_required(VERSION 3.14)
project(SabreRecon)

find_package(ROOT REQUIRED)
set(SABRERECON_BINARY_DIR ${CMAKE_CURRENT_SOURCE_DIR}/bin)
set(SABRERECON_LIBRARY_DIR ${CMAKE_CURRENT_SOURCE_DIR}/lib)

set(CMAKE_CXX_STANDARD 17)

add_subdirectory(src/vendor/catima)
add_subdirectory(src)
