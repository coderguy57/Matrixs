cmake_minimum_required(VERSION 3.0.0)

add_executable(tests_project "${PROJECT_SOURCE_DIR}/tests/test_main.cpp")

include(FetchContent)
FetchContent_Declare(
        Catch2
        GIT_REPOSITORY https://github.com/catchorg/Catch2
        GIT_TAG v2.13.10
)

FetchContent_MakeAvailable(Catch2)

target_link_libraries(tests_project
    PRIVATE
    Catch2::Catch2
)

# Add memory leak detection to all tests
add_test(NAME MemoryLeaks COMMAND valgrind --leak-check=full --error-exitcode=1 $<TARGET_FILE:tests_project>)

list(APPEND CMAKE_MODULE_PATH ${catch2_SOURCE_DIR}/contrib)

target_include_directories(tests_project PUBLIC
	${PROJECT_SOURCE_DIR}/src)

include(Catch)
catch_discover_tests(tests_project)