cmake_minimum_required(VERSION 3.15)
project(MyProject)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_FLAGS "-O3")
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# Include directories
include_directories(${PROJECT_SOURCE_DIR}/include)

# Add an option to build Python bindings
option(BUILD_PYTHON_BINDINGS "Build Python bindings" OFF)

# Add the executable for the main program
add_executable(main_program src/main.cpp)

# Add the executable for the LU decomposition example
add_executable(lu_example examples/lu_example.cpp)

# Define the path to the test directory for the LU example
target_compile_definitions(lu_example PRIVATE TEST_DIR="${PROJECT_SOURCE_DIR}/test")

# Link to the Accelerate framework on macOS
if(APPLE)
    target_link_libraries(main_program "-framework Accelerate")
    target_link_libraries(lu_example "-framework Accelerate")
endif()

# Conditionally add pybind11 and the Python module
if(BUILD_PYTHON_BINDINGS)
    # Use the new FindPython3 module
    find_package(Python3 REQUIRED COMPONENTS Interpreter Development)
    find_package(pybind11 REQUIRED)

    include_directories(${Python3_INCLUDE_DIRS})
    pybind11_add_module(DNR src/bindings.cpp)
    target_link_libraries(DNR PRIVATE ${Python3_LIBRARIES})
endif()

# Print the flags to verify
message(STATUS "CMake CXX Flags: ${CMAKE_CXX_FLAGS}")
message(STATUS "Build Python bindings: ${BUILD_PYTHON_BINDINGS}")

