#ifndef OPENFLIPPER_SMOOTHING_HH
#define OPENFLIPPER_SMOOTHING_HH

#include <ObjectTypes/TriangleMesh/TriangleMesh.hh>

class Smoothing
{
private:
    Smoothing();
public:
    static void smooth(TriMesh& _mesh, OpenMesh::SmartVertexHandle _vh);

};

#endif // OPENFLIPPER_SMOOTHING_HH
