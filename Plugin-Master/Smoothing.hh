#ifndef OPENFLIPPER_SMOOTHING_HH
#define OPENFLIPPER_SMOOTHING_HH

#include <ObjectTypes/TriangleMesh/TriangleMesh.hh>

class Smoothing
{
public:
    Smoothing(TriMesh& _mesh):mesh_(_mesh){}
public:
    void smooth(OpenMesh::SmartVertexHandle _vh);
private:
    void phi();
private:
    TriMesh& mesh_;
};

#endif // OPENFLIPPER_SMOOTHING_HH
