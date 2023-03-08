#ifndef OPENFLIPPER_VERTEXDISPLACEMENT_HH
#define OPENFLIPPER_VERTEXDISPLACEMENT_HH

#include <ObjectTypes/TriangleMesh/TriangleMesh.hh>

class VertexDisplacement {
public:
    VertexDisplacement(TriMesh& _mesh): mesh_(_mesh){}
    ~VertexDisplacement(){}
public:
    void displace(const ACG::Vec3d _displacement, const int _constraint_vh);
private:
    TriMesh& mesh_;

};
#endif // VERTEXDISPLACEMENT_HH
