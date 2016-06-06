cmake_minimum_required(VERSION 3.5)
project(EMOCT)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
set(EXEC_NAME EMOCT)
SET(QTCP_PROJECT_DIR .)

set(SOURCE_FILES
        gact_main.cpp
        CUDAMath.h
        MersenneTwister.h
        GACT.cpp
        mtrand.cpp
        ilab.h
        ilab/art.h
        ilab/distribution.h
        ilab/projection.h
        ilab/projector.h
         )
add_executable(EMOCT ${SOURCE_FILES})