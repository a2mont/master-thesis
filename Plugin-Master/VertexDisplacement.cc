#include "VertexDisplacement.hh"

void VertexDisplacement::displace(const ACG::Vec3d _displacement, const int _constraint_vh){

    std::cout << "Displacement: " << _displacement << std::endl;
    for(auto vh: mesh_.vertices()){
        if (vh.idx() == _constraint_vh) {
            std::cout << "Constraint vertex id: " << vh.idx() << std::endl;
            auto pt = mesh_.point(vh);
            pt += _displacement;
            mesh_.point(vh) = pt;
        }
    }

}
