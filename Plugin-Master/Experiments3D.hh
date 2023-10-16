#ifndef OPENFLIPPER_EXPERIMENTS3D_HH
#define OPENFLIPPER_EXPERIMENTS3D_HH

#include "TetLoop.hh"
#include "TetrahedralizedVoxelGridGenerator.hh"

class Experiment3D
{
    using Point = ACG::Vec3d;
public:
    Experiment3D(TetrahedralMesh& _mesh, const double _q_min, std::map<int,int>& _constraint_vhs):
        mesh_(_mesh),
        loop_(_mesh, _q_min, _constraint_vhs, true)
    {}
private:
    TetrahedralMesh& mesh_;
    TetLoop loop_;

public:
    //3D
    void generate_torsion_mesh(double torsion_turns_count);

};

#endif // EXPERIMENTS_HH
