# Specify the minimum CMAKE version required
CMAKE_MINIMUM_REQUIRED(VERSION 2.8)

# Your project's name
PROJECT(superradiance)

# Instruct CMake to inspect the following subfolders
if(UNIX)
    SET(CMAKE_CXX_COMPILER i686-w64-mingw32-gcc)
elseif(WIN32)
    INCLUDE_DIRECTORIES("C:/MinGW/include")
endif()

SET(CMAKE_CXX_FLAGS "-std=c++11 -Wimplicit-function-declaration")

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

SET(SWSH_SOURCE swsh/swsh.c swsh/matrix3d.c swsh/matrix3d.h)
ADD_EXECUTABLE(swsh ${SWSH_SOURCE})
TARGET_LINK_LIBRARIES(swsh m gsl gslcblas)