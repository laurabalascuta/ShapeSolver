ShapeSolver.exe: ShapeSolver.cpp Entities\Point.cpp Entities\Vector.cpp Entities\Triangle.cpp Entities\TriangleFactory.cpp Solver\OrangeTriangleSolver.cpp Utils\IO.cpp Utils\Math.cpp
	cl /nologo /EHsc $** 

clean: dummy
	-@del ShapeSolver.exe
	-@del ShapeSolver.obj

dummy:
