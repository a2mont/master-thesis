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
    static bool t_chebyshev_centroid();
    static bool t_custom_chebyshev_centroid(std::string _filename="mesh_dump0_3D.ovm");
private:
    Tests();
    inline static const std::string LOGS_MESH = "../../../../Plugin-Master/logs/meshes/";
};

#endif // TESTS_HH
