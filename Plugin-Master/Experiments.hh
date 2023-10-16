#ifndef OPENFLIPPER_EXPERIMENTS_HH
#define OPENFLIPPER_EXPERIMENTS_HH

#include <ObjectTypes/TriangleMesh/TriangleMesh.hh>
#include "TetLoop.hh"
#include "TriangleLoop.hh"

class Experiment
{
    using Point = ACG::Vec3d;
public:
    Experiment(TriMesh& _mesh, TriMesh& _worldMesh,std::map<int,Point> _basePoints, const double _q_min, std::map<int,int>& _constraint_vhs):
        mesh_(_mesh),
        loop(_mesh, _worldMesh, _q_min, _constraint_vhs, true),
        basePoints_(_basePoints)
    {}

private:
    TriMesh& mesh_;
    TriangleLoop loop;
    std::map<int,Point> basePoints_;
    std::map<int,std::vector<Point>> targetPoints_;
    bool initialized_ = false;
    std::vector<OpenMesh::SmartVertexHandle> leftBoundary_;
    std::vector<OpenMesh::SmartVertexHandle> rightBoundary_;
    std::vector<OpenMesh::SmartVertexHandle> topBoundary_;
    std::vector<OpenMesh::SmartVertexHandle> bottomBoundary_;


    void updateBoundaries();
    Eigen::Vector3f generateParabola(Point _a, Point _b, Point _c);
    Point findMidPoint(Point _position, std::vector<OpenMesh::SmartVertexHandle> _boundary);
public:
    //2D
    void bend2D(const int _timesteps, const int _currentT);
    void stretch2D(const int _timesteps, const int _currentT);
    void compress2D(const int _timesteps, const int _currentT);
    //3D
    void generate_torsion_mesh(double torsion_turns_count);

};

#endif // EXPERIMENTS_HH
