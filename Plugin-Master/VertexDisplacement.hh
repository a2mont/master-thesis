#ifndef OPENFLIPPER_VERTEXDISPLACEMENT_HH
#define OPENFLIPPER_VERTEXDISPLACEMENT_HH

#include <ObjectTypes/TriangleMesh/TriangleMesh.hh>

class VertexDisplacement {
private:
    VertexDisplacement();
public:
    static void displace(TriMesh& _mesh, const ACG::Vec3d _displacement, const int _constraint_vh, const bool _verbose = false);

};
#endif // OPENFLIPPER_VERTEXDISPLACEMENT_HH
