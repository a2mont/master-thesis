#ifndef TESTS_HH
#define TESTS_HH

#include "TetLoop.hh"


class Tests
{

public:
    static bool t_EdgeRemoval();
    static bool t_FaceRemoval();
    static bool t_EdgeContraction();
    static bool runAll();
private:
    Tests();
    static void computeQuality(TetLoop::PriorityQueue &_queue, TetrahedralMesh &_mesh);
};

#endif // TESTS_HH
