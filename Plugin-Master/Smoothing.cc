#include "Smoothing.hh"

void Smoothing::smooth(TriMesh& _mesh, OpenMesh::SmartVertexHandle _vh){
    ACG::Vec3d pt(0.,0.,0.);
    unsigned int n_neighbours = 0;

    if(_vh.is_boundary())
        return;
    for(auto vvh: _vh.vertices()){
        pt += _mesh.point(vvh);
        n_neighbours++;
    }
    pt /= n_neighbours;
    _mesh.point(_vh) = pt;
}

void Smoothing::smooth(TetrahedralMesh& _mesh, OpenVolumeMesh::VertexHandle _vh){
    ACG::Vec3d pt(0.,0.,0.);
    unsigned int n_neighbours = 0;

    if(_mesh.is_boundary(_vh))
        return;
    for(auto vvh: _mesh.vertex_vertices(_vh)){
        pt += _mesh.vertex(vvh);
        n_neighbours++;
    }
    pt /= n_neighbours;
    _mesh.set_vertex(_vh, pt);
}

