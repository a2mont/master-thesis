#ifndef OPENFLIPPER_SMOOTHING_HH
#define OPENFLIPPER_SMOOTHING_HH

#include <ObjectTypes/TriangleMesh/TriangleMesh.hh>
#include <ObjectTypes/TetrahedralMesh/TetrahedralMesh.hh>

class Smoothing
{
private:
    Smoothing();
public:
    static void smooth(TriMesh& _mesh, OpenMesh::SmartVertexHandle _vh);
    static void smooth(TetrahedralMesh& _mesh, OpenVolumeMesh::VertexHandle _vh);

};

#endif // OPENFLIPPER_SMOOTHING_HH
