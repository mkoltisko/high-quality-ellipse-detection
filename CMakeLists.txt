cmake_minimum_required(VERSION 3.17)
project(Help_David)

set(CMAKE_CXX_STANDARD 14)

include_directories(.)

add_executable(Help_David
        Cluster.h
        generateEllipseCandidates.cpp
        Gradient.h
        HelperFunctions.h
        LineSegment.h
        MagicNumbers.h
        NFA.h
        Region.h
        Structs.h)
