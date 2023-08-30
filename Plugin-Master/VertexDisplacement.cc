#include "VertexDisplacement.hh"
using namespace OpenVolumeMesh;

void VertexDisplacement::displace(TriMesh& _mesh, const ACG::Vec3d _displacement,
                                  const std::map<int,int> _constraint_vhs, const bool _verbose){
    if(_verbose)
        std::cout << "Displacement: " << _displacement << std::endl;

    for(auto vh_id: _constraint_vhs){
        auto vh = _mesh.vertex_handle(vh_id.first);
        if(_verbose)
            std::cout << "Constraint vertex id: " << vh.idx() << std::endl;
        auto pt = _mesh.point(vh);
        pt += _displacement;
        _mesh.point(vh) = pt;
    }

}
void VertexDisplacement::displace(TetrahedralMesh& _mesh, const ACG::Vec3d _displacement,
                                  const std::map<int,int> _constraint_vhs, const bool _verbose){
    if(_verbose)
        std::cout << "Displacement: " << _displacement << std::endl;

    for(auto vh_id: _constraint_vhs){
        auto vh = VertexHandle(vh_id.first);
        if(_verbose)
            std::cout << "Constraint vertex id: " << vh.idx() << std::endl;
        auto pt = _mesh.vertex(vh);
        pt += _displacement;
        _mesh.set_vertex(vh, pt);
    }

}
