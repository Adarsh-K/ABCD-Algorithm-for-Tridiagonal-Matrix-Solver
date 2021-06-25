cmake_minimum_required(VERSION 3.14...3.19)

project(abcd LANGUAGES CXX)

find_package(Armadillo REQUIRED)

include_directories(${ARMADILLO_INCLUDE_DIRS})

set(SOURCES
		src/main.cpp
		src/partition.cpp
		src/read_matrix.cpp
		src/get_matrices.cpp
		src/utils.cpp
		)

add_executable(${PROJECT_NAME} ${SOURCES})

target_link_libraries(${PROJECT_NAME} ${ARMADILLO_LIBRARIES})
