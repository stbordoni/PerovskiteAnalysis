cmake_minimum_required(VERSION 3.10)
project(PerovskiteAna)

set(CMAKE_CXX_STANDARD 11)

option( WITH_ROOT "Build also binaries depending on ROOT" ON )

include( ${CMAKE_SOURCE_DIR}/cmake/cmessage.cmake )
include( ${CMAKE_SOURCE_DIR}/cmake/dependencies.cmake )




find_package(ROOT CONFIG REQUIRED Graf Gpad Hist RIO Tree)

add_subdirectory( ${CMAKE_SOURCE_DIR}/src)