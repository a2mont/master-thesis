#ifndef OPENFLIPPER_SMOOTHING_HH
#define OPENFLIPPER_SMOOTHING_HH

#include <ObjectTypes/TriangleMesh/TriangleMesh.hh>
#include <ObjectTypes/TetrahedralMesh/TetrahedralMesh.hh>
#include "ortools/linear_solver/linear_solver.h"

class Smoothing
{
private:
    Smoothing();
    static void triangle_normal_and_centroid(const TetrahedralMesh &_mesh, const OpenVolumeMesh::HalfFaceHandle &hf, ACG::Vec3d &normal, ACG::Vec3d &centroid, bool printDebug=false);
    static bool find_chebyshev_center(const TetrahedralMesh &mesh, const std::set<OpenVolumeMesh::HalfFaceHandle> &constraint_hfs, const double radius_lower_bound, ACG::Vec3d &new_position, bool printDebug=false);
public:
    static void smooth(TriMesh& _mesh, OpenMesh::SmartVertexHandle _vh);
    static void smooth(TetrahedralMesh& _mesh, OpenVolumeMesh::VertexHandle _vh);

};

#endif // OPENFLIPPER_SMOOTHING_HH
