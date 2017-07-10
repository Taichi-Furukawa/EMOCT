cmake_minimum_required(VERSION 3.5)
project(EMOCT)

set(EXEC_NAME EMOCT)
SET(QTCP_PROJECT_DIR .)
set(CMAKE_CXX_COMPILER /usr/local/opt/llvm/bin/clang++)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -fopenmp -L/usr/local/opt/llvm/lib -Wl,-rpath,/usr/local/opt/llvm/lib")

set(SOURCE_FILES
        gact_main.cpp
        CUDAMath.h
        MersenneTwister.h
        mtrand.cpp
        ilab.h
        ilab/art.h
        ilab/distribution.h
        ilab/projection.h
        ilab/projector.h
        EMO.cpp EMO.h
        InverseDomain.cpp InverseDomain.h)

include_directories(src)

set(EIGEN_INCLUDE_PATH /usr/local/include/eigen3/)
set(EIGEN_LIBRARY_PATH /usr/local/include/eigen3/Eigen)

set(BOOST_INCLUDE_PATH /usr/local/Cellar/boost/1.64.0_1/include)
set(BOOST_LIBRARY_PATH /usr/local/Cellar/boost/1.64.0_1/lib)
find_package(OpenMP)
if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()

find_package(Eigen3 COMPONENTS system filesystem REQUIRED)
find_package(BOOST COMPONENTS system filesystem REQUIRED)
include_directories(${EIGEN_INCLUDE_PATH})
include_directories(${BOOST_INCLUDE_PATH})

add_executable(EMOCT ${SOURCE_FILES})
find_package(OpenCV REQUIRED)
include_directories(${OpenCV_INCLUDE_DIRS})

target_link_libraries(EMOCT ${OpenCV_LIBS})