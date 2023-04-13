#include "Highlight.hh"

void Highlight::highlight_triangles(TriMesh& _mesh, const bool _verbose){
    for(auto fh: _mesh.faces()){
        double quality = QualityEvaluation::evaluate(fh, _mesh, _verbose);
        _mesh.set_color(fh, ACG::Vec4f(1-quality, quality, 0, 1));
    }

}
