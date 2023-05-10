#include "VertexDisplacement.hh"

void VertexDisplacement::displace(TriMesh& _mesh, const ACG::Vec3d _displacement, const int _constraint_vh, const bool _verbose){
    if(_verbose)
        std::cout << "Displacement: " << _displacement << std::endl;
    for(auto vh: _mesh.vertices()){
        if (vh.idx() == _constraint_vh) {
            if(_verbose)
                std::cout << "Constraint vertex id: " << vh.idx() << std::endl;
            auto pt = _mesh.point(vh);
            pt += _displacement;
            _mesh.point(vh) = pt;
        }
    }

}
