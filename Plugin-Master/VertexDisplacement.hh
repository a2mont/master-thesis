#ifndef OPENFLIPPER_VERTEXDISPLACEMENT_HH
#define OPENFLIPPER_VERTEXDISPLACEMENT_HH

#include <ObjectTypes/TriangleMesh/TriangleMesh.hh>
#include <ObjectTypes/TetrahedralMesh/TetrahedralMesh.hh>

class VertexDisplacement {
private:
    VertexDisplacement();
public:
    static void displace(TriMesh& _mesh, const ACG::Vec3d _displacement,
                         const std::map<int,int> _constraint_vhs, const bool _verbose = false);
    static void displace(TetrahedralMesh& _mesh, const ACG::Vec3d _displacement,
                         const std::map<int,int> _constraint_vhs, const bool _verbose = false);

};
#endif // OPENFLIPPER_VERTEXDISPLACEMENT_HH
