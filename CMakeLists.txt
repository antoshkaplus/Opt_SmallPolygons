# cmake -DCMAKE_BUILD_TYPE=Debug
# cmake -DCMAKE_BUILD_TYPE=Release

cmake_minimum_required (VERSION 2.8.11)
project (SmallPolygon)

list(APPEND CMAKE_CXX_FLAGS "-std=c++1y")

set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -d -Wall")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3 -Wall")

set(BinDir "../bin")
set(LibDir "../lib")

set( CMAKE_RUNTIME_OUTPUT_DIRECTORY ${BinDir} )
set( CMAKE_LIBRARY_OUTPUT_DIRECTORY ${LibDir} )
set( CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${LibDir} )
foreach( OUTPUTCONFIG ${CMAKE_CONFIGURATION_TYPES} )
    string( TOUPPER ${OUTPUTCONFIG} OUTPUTCONFIG )
    set( CMAKE_RUNTIME_OUTPUT_DIRECTORY_${OUTPUTCONFIG} ${BinDir} )
    set( CMAKE_LIBRARY_OUTPUT_DIRECTORY_${OUTPUTCONFIG} ${LibDir} )
    set( CMAKE_ARCHIVE_OUTPUT_DIRECTORY_${OUTPUTCONFIG} ${LibDir} )
endforeach( OUTPUTCONFIG CMAKE_CONFIGURATION_TYPES )

file(GLOB SOURCE_FILES "./src/*.cpp")
file(GLOB HEADER_FILES "./include/*.hpp")

file(GLOB SUBMIT_FILES "./submit/version.cpp")

include_directories("./include")
include_directories("/usr/local/include")
include_directories("/usr/include")

find_library(ANT_LIBRARY ant)
if (NOT ANT_LIBRARY) 
    message(FATAL_ERROR "ant library not found")    
endif()

add_library(small_polygons ${SOURCE_FILES} ${HEADER_FILES})
target_link_libraries(small_polygons ${ANT_LIBRARY})
    
add_executable(solver "./app/solver.cpp")
target_link_libraries(solver small_polygons)

add_executable(test "./app/test.cpp")
target_link_libraries(test small_polygons ${ANT_LIBRARY})

add_library(submit STATIC ${SUBMIT_FILES})

