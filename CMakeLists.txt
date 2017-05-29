# Specify the minimum CMAKE version required
CMAKE_MINIMUM_REQUIRED(VERSION 3.0.0)

PROJECT(superradiance VERSION 0.0.2)


# Set output directories
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/../lib")
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/../lib")
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/../bin")


# OS specific variables
if(UNIX)
    SET(CMAKE_CXX_FLAGS "-std=c++11 -Wimplicit-function-declaration")
    SET(OS_INCLUDE_DIRS "/usr/include" "/usr/local/include")
    SET(OS_LIB_DIRS "/usr/lib" "/usr/local/lib")
    #SET(CMAKE_C_COMPILER "x86_64-w64-mingw32-gcc")
elseif(WIN32)
    SET(CMAKE_CXX_FLAGS "-posix -std=c++11 -Wimplicit-function-declaration")
    SET(OS_INCLUDE_DIRS "C:/MinGW/include")    
endif()

INCLUDE_DIRECTORIES("${CMAKE_BINARY_DIR}/../include" ${OS_INCLUDE_DIRS})
LINK_DIRECTORIES("${CMAKE_BINARY_DIR}/../lib" ${OS_LIB_DIRS})


# Source related variables
SET(TENSOR3_SRC
    src/gsl/gsl_tensor3.c
    src/gsl/gsl_tensor3_long_double.c           
)
SET(SWSH_SRC
    src/test.c
    src/diffeq.c
    src/swsh.c
)


# First build the libraries
ADD_LIBRARY(gsl-tensor3 STATIC ${TENSOR3_SRC})


# Sencodly build/link the execuatables (order is important)
ADD_EXECUTABLE(swsh ${SWSH_SRC})
TARGET_LINK_LIBRARIES(swsh gsl-tensor3 m gsl gslcblas)