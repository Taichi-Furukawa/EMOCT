cmake_minimum_required(VERSION 3.5)
project(EMOCT)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
set(EXEC_NAME EMOCT)
SET(QTCP_PROJECT_DIR .)

set(SOURCE_FILES
        ART_Main.cpp
        ProjectionData.cpp
        ART.cpp
        CTUtility.cpp
        DensityDistribution.cpp
        CUDAMath.h
        MersenneTwister.h
        Xorshift.h
        )
add_executable(EMOCT ${SOURCE_FILES})