#include "Highlight.hh"

ACG::Vec4f colors [4] = {ACG::Vec4f(1,0,0,1), ACG::Vec4f(1,0.5,0,1), ACG::Vec4f(0,1,0,1), ACG::Vec4f(0,0,1,1)};

void Highlight::highlight_worst_triangles(TriMesh& _mesh, double _treshold){
    QualityEvaluation qe(_mesh);
    std::vector<QualityEvaluation::Triangle> worstTriangles;
    for(auto fh: _mesh.faces()){
        qe.evaluate(fh);
    }
    auto queue = qe.get_face_quality_queue();
    auto worstTriangle = queue.top();
    int worstTriId = worstTriangle.face_id_;
    while (worstTriangle.quality_ < _treshold) {
        worstTriangles.emplace_back(worstTriangle);
        queue.pop();
        worstTriangle = queue.top();
    }

    for(auto tri: worstTriangles){
        _mesh.set_color(_mesh.face_handle(tri.face_id_), colors[1]);
    }
    _mesh.set_color(_mesh.face_handle(worstTriId), colors[0]);

    qe.reset_queue();
}
