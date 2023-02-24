cmake_minimum_required(VERSION 3.0.0)
project(Matrixs VERSION 0.1.0)

include(CTest)
enable_testing()

add_library(Matrixs src/matrix.cpp)

set(CPACK_PROJECT_NAME ${PROJECT_NAME})
set(CPACK_PROJECT_VERSION ${PROJECT_VERSION})
include(CPack)


if (BUILD_TESTING)
    enable_testing()
	add_subdirectory(tests)
endif()
