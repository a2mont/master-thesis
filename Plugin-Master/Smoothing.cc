#include "Smoothing.hh"

void Smoothing::smooth(OpenMesh::SmartVertexHandle _vh){
    ACG::Vec3d pt(0.,0.,0.);
    unsigned int n_neighbours = 0;

    if(_vh.is_boundary())
        return;
    for(auto vvh: _vh.vertices()){
        pt += mesh_.point(vvh);
        n_neighbours++;
    }
    pt /= n_neighbours;
    mesh_.point(_vh) = pt;
    // compute phi_0 and A_0, phi cf below, A_0 the initial active set
}


void Smoothing::phi(){
    // min theta_j, j in S_i the set of all angles in submesh,
}
