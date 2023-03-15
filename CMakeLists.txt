cmake_minimum_required(VERSION 3.0.0)

project(HdivStokes)

include_directories(
	"${PROJECT_SOURCE_DIR}/mhstokes"
)

find_package(NeoPZ REQUIRED HINTS ${CMAKE_SOURCE_DIR}/../neopz_install/ ${CMAKE_SOURCE_DIR}/neopz_install/)

file(GLOB_RECURSE PZ_SRC
     "../neopz/*.h"
     "../neopz/*.cpp"
)

add_custom_target( GlobalFiles SOURCES ${PZ_SRC})

add_subdirectory(mhstokes)
add_subdirectory(lib)
add_subdirectory(executable)
