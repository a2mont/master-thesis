#ifndef TESTS_HH
#define TESTS_HH

#include "TetLoop.hh"
#include "TetrahedralizedVoxelGridGenerator.hh"
#include "Experiments3D.hh"

using namespace OpenVolumeMesh;
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
    static bool t_flip23();
    static bool t_custom_EdgeRemoval(std::string _filename);
    static bool t_flip32();
    static bool t_multiface();
    static bool t_custom_multiface(std::string _filename="dump_multiface.ovm");
    static bool t_quality_evaluation();
    static bool t_speed();
private:
    Tests();
    inline static const std::string LOGS_MESH = "../../../../Plugin-Master/logs/meshes/";
};

#endif // TESTS_HH
