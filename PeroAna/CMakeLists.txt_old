cmake_minimum_required(VERSION 3.10)
project(PerovskiteAna)

set(CMAKE_CXX_STANDARD 11)

# Add the path to ROOT libraries and include directories
list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS})
find_package(ROOT REQUIRED COMPONENTS Core RIO Tree)

include_directories(${ROOT_INCLUDE_DIRS})
add_definitions(${ROOT_DEFINITIONS})
list(APPEND EXTRA_LIBS ${ROOT_LIBRARIES})


# Set the source directory
set(SOURCE_DIR src)

# Add the source file in the src directory
#set(SOURCES ${SOURCE_DIR}/readData.cpp)
#add_executable(readData ${SOURCES})


#Add the executable for readData.cpp
add_executable(readData ${SOURCE_DIR}/readData.cpp)


# Add executable for MyAna
add_executable(PerovAna ${SOURCE_DIR}/PerovAna.cpp ${SOURCE_DIR}/event.cpp)

# Set the output directory for the executable
#set_target_properties(readData PROPERTIES
#    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/build"
#)

# Set the output directory for the executable
#set_target_properties(PerovAna PROPERTIES
#    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/build"
#)

#target_link_libraries(readData ${ROOT_LIBRARIES})
#target_link_libraries(PerovAna ${ROOT_LIBRARIES})
target_link_libraries(readData ${EXTRA_LIBS})
target_link_libraries(PerovAna ${EXTRA_LIBS})