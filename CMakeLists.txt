cmake_minimum_required(VERSION 3.15)
project(untitled)

set(CMAKE_CXX_STANDARD 17)
SET(CMAKE_BUILD_TYPE Debug)
set(CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake/modules)
message(STATUS "CMAKE_MODULE_PATH: ${CMAKE_MODULE_PATH}")
add_subdirectory(include/fmt EXCLUDE_FROM_ALL)
# add_subdirectory(include/vorbis EXCLUDE_FROM_ALL)

find_package(PNG REQUIRED)
# Find libsndfile
find_library(SNDFILE_LIBRARIES NAMES sndfile)

# Check if libsndfile was found
if(NOT SNDFILE_LIBRARIES)
    message(FATAL_ERROR "Could not find libsndfile")
endif()

if (VORBIS_LIBRARIES)
    message(STATUS "Found libvorbis: ${VORBIS_LIBRARIES}")
endif ()


add_executable(untitled main.cpp)
target_link_libraries(untitled ${PNG_LIBRARY} ${SNDFILE_LIBRARIES} fmt::fmt-header-only)

target_include_directories(untitled PRIVATE
        ${PROJECT_SOURCE_DIR}/include/fmt/include
        ${PROJECT_SOURCE_DIR}/include/libsndfile
        ${PNG_INCLUDE_DIR}
        )
