cmake_minimum_required(VERSION 3.17)
project(ShapeSolver)

set(CMAKE_CXX_STANDARD 17)

add_executable(ShapeSolver main.cpp Entities/Entity.h
        Entities/Point.cpp Entities/Point.h
        Entities/Triangle.cpp Entities/Triangle.h
        Entities/Vector.cpp Entities/Vector.h
        Solver/ShapeSolver.h
        Solver/OrangeTriangleSolver.cpp Solver/OrangeTriangleSolver.h
        Utils/IO.cpp Utils/IO.h
        Utils/Math.cpp Utils/Math.h Entities/Factory.h Entities/TriangleFactory.cpp Entities/TriangleFactory.h)