cmake_minimum_required(VERSION 3.15)
project(untitled)

set(CMAKE_CXX_STANDARD 17)
SET(CMAKE_BUILD_TYPE Debug)
set(CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake/modules)
message(STATUS "CMAKE_MODULE_PATH: ${CMAKE_MODULE_PATH}")
add_subdirectory(include/fmt EXCLUDE_FROM_ALL)
# add_subdirectory(include/vorbis EXCLUDE_FROM_ALL)

find_package(PNG REQUIRED)

set(VORBISFILE_LIB "/usr/local/Cellar/libvorbis/1.3.7/lib/libvorbisfile.a")
set(VORBISENC_LIB "/usr/local/Cellar/libvorbis/1.3.7/lib/libvorbisenc.a")
set(VORBIS_LIB "/usr/local/Cellar/libvorbis/1.3.7/lib/libvorbis.a")
set(OGG_LIB "/usr/local/Cellar/libogg/1.3.5/lib/libogg.a")
# find_package(Vorbis REQUIRED)

if (VORBIS_LIBRARIES)
    message(STATUS "Found libvorbis: ${VORBIS_LIBRARIES}")
endif ()


add_executable(untitled main.cpp)
target_link_libraries(untitled ${PNG_LIBRARY} ${OGG_LIB} ${VORBIS_ENC} ${VORBIS_LIB} ${VORBISFILE_LIB} fmt::fmt-header-only)

target_include_directories(untitled PRIVATE
        ${PROJECT_SOURCE_DIR}/include/fmt/include
        ${PROJECT_SOURCE_DIR}/include/vorbis/include
        "/usr/local/Cellar/libogg/1.3.5/include"
        ${PNG_INCLUDE_DIR}
        )
