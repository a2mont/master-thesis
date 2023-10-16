#ifndef TESTS_HH
#define TESTS_HH

#include "TetLoop.hh"
#include "TetrahedralizedVoxelGridGenerator.hh"


class Tests
{

public:
    static bool t_EdgeRemoval();
    static bool t_FaceRemoval();
    static bool t_EdgeContraction();
    static bool t_StressEdgeRemoval();
    static bool runAll();
    static bool t_StressFaceRemoval();
private:
    Tests();
};

#endif // TESTS_HH
