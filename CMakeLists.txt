cmake_minimum_required(VERSION 3.16.3)

# Projects
project(GLOCR VERSION 0.0.1)
set(PHASE1 phase1_tests)
set(PHASE2 phase2_tests)
set(IMGPROC imageprocessing_tests)
set(UTILS utils_tests)
set(IAMDOCCLEANER cleanIAMDocs_tool)

# Standard 20
set(CMAKE_CXX_STANDARD 20)
message(STATUS "CXX Standard: ${CMAKE_CXX_STANDARD}")

# Project Structure Info
include_directories(./include)
link_directories(./libs)
add_subdirectory(src)

# Move data dependences to the binary directory
configure_file(
    ${CMAKE_SOURCE_DIR}/data/document.png 
    ${CMAKE_CURRENT_BINARY_DIR}/data/document.png
    COPYONLY
)
configure_file(
    ${CMAKE_SOURCE_DIR}/data/document_cleaned.png 
    ${CMAKE_CURRENT_BINARY_DIR}/data/document_cleaned.png
    COPYONLY
)
configure_file(
    ${CMAKE_SOURCE_DIR}/data/word.png 
    ${CMAKE_CURRENT_BINARY_DIR}/data/word.png
    COPYONLY
)

# OpenMP
find_package(OpenMP)
if (OPENMP_FOUND)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
    
    set(OMP_COMPILE_OPTS ${OpenMP_CXX_FLAGS})
    set(OMP_LINK_OPTS OpenMP::OpenMP_CXX)

    message(WARNING "OpenMP found, but the multithreading may not be enabled. This is a bug on my part. Please check README.md")
endif() 

# Phase 1
add_executable(${PHASE1} tests/phase1_test.cpp)
target_compile_options(${PHASE1} PRIVATE -Wall ${OMP_COMPILE_OPTS})
target_link_libraries(${PHASE1} PUBLIC utils imageprocessing imagesegmentation png z ${OMP_LINK_OPTS})

# Phase 2
add_executable(${PHASE2} tests/phase2_test.cpp)
target_compile_options(${PHASE2} PRIVATE -Wall ${OMP_COMPILE_OPTS})
target_link_libraries(${PHASE2} PUBLIC utils imageprocessing imagesegmentation png z ${OMP_LINK_OPTS})

# Image Processing
add_executable(${IMGPROC} tests/imageprocessing_test.cpp)
target_compile_options(${IMGPROC} PRIVATE -Wall ${OMP_COMPILE_OPTS})
target_link_libraries(${IMGPROC} PUBLIC utils imageprocessing png z ${OMP_LINK_OPTS})

# IAM Document Cleaning Tool
add_executable(${IAMDOCCLEANER} src/cleanIAMDocs.cpp)
target_compile_options(${IAMDOCCLEANER} PRIVATE -Wall ${OMP_COMPILE_OPTS})
target_link_libraries(${IAMDOCCLEANER} PUBLIC utils imageprocessing png z ${OMP_LINK_OPTS})

# Utils
add_executable(${UTILS} tests/utils_test.cpp)
target_link_libraries(${UTILS} utils)