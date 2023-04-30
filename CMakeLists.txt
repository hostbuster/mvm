cmake_minimum_required(VERSION 3.15)
project(untitled)

set(CMAKE_CXX_STANDARD 17)
SET(CMAKE_BUILD_TYPE Debug)
set(CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake/modules)
message(STATUS "CMAKE_MODULE_PATH: ${CMAKE_MODULE_PATH}")
add_subdirectory(include/fmt EXCLUDE_FROM_ALL)
add_subdirectory(include/Gist EXCLUDE_FROM_ALL)
# add_subdirectory(include/vorbis EXCLUDE_FROM_ALL)

find_package(PNG REQUIRED)
# Find libsndfile
find_library(SNDFILE_LIBRARIES NAMES sndfile)

add_library(kiss_fft STATIC include/Gist/libs/kiss_fft130/kiss_fft.c include/Gist/libs/kiss_fft130/kiss_fft.h)
# target_include_directories(kiss_fft PUBLIC libs/kiss_fft)


# Check if libsndfile was found
if(NOT SNDFILE_LIBRARIES)
    message(FATAL_ERROR "Could not find libsndfile")
endif()


add_executable(untitled main.cpp)
target_link_libraries(untitled ${PNG_LIBRARY} ${SNDFILE_LIBRARIES} fmt::fmt-header-only Gist kiss_fft)

target_include_directories(untitled PRIVATE
        ${PROJECT_SOURCE_DIR}/include/fmt/include
        ${PROJECT_SOURCE_DIR}/include/libsndfile
        ${PROJECT_SOURCE_DIR}/include/Gist/src
        ${PROJECT_SOURCE_DIR}/include/Gist/libs/kiss_fft130
        ${PNG_INCLUDE_DIR}
        )
