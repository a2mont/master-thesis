#include "MainLoop.hh"

void MainLoop::loop(ACG::Vec3d _displacement, int _constraint_vh, bool _verbose,int _max_iter ){

    for(int i=0; i < _max_iter; i++){
        QualityEvaluation::PriorityQueue worstTriangles;
        QualityEvaluation::Triangle worstTriangle(0,0);

        vd_.displace(_displacement, _constraint_vh, _verbose);

        for(auto fh: mesh_.faces()){
            qe_.evaluate(fh, _verbose);
        }

        worstTriangle = qe_.get_face_quality_queue().top();
        while(worstTriangle.quality_ < q_min){
            worstTriangles.push(worstTriangle);
            qe_.get_face_quality_queue().pop();
            worstTriangle = qe_.get_face_quality_queue().top();
        }

        if(worstTriangles.empty())
            std::cout << "All triangles are of good enough quality\nWorst triangle: " << worstTriangle.toString() << std::endl;
        else
            std::cout << worstTriangles.size() << " triangles of bad quality"
                      << "\nWorst triangle: " << worstTriangles.top().toString()
                      << std::endl;
        improve_mesh(worstTriangles);

        qe_.reset_queue();
    }
}

void MainLoop::improve_mesh(QualityEvaluation::PriorityQueue _badTriangles){
    // B <- set of triangles in M with q < q_min
    // foreach t in B
    // if t still exists and q < q_min improve_tet(M,t,q_min)
    while (!_badTriangles.empty()) {
        auto t = _badTriangles.top();
        _badTriangles.pop();
        // TODO check if t exists mesh_isdeleted)
        if(!mesh_.status(mesh_.face_handle(t.face_id_)).deleted())
            improve_triangle(t);
    }
}

void MainLoop::improve_triangle(QualityEvaluation::Triangle _t){
    // A <- t, a triangle from the mesh
    QualityEvaluation::PriorityQueue A;
    A.push(_t);
    for (int i = 0; i < 10; i++) {
        bool changed = false;
        do {
            // A <- topological_pass(A,M)
            changed = topologial_pass(&A);
            if(A.top().quality_ >= q_min)
                return;
        } while (changed);

        // A <- edge_contraction_pass(A,M)
        edge_contraction_pass(&A);
        if(A.top().quality_ >= q_min)
            return;

        // A <- insertion_pass(A,M)
        insertion_pass(&A);
        if(A.top().quality_ >= q_min)
            return;

        // smoothing_pass()
        smoothing_pass(&A);
        if(A.top().quality_ >= q_min)
            return;
    }
}

bool MainLoop::topologial_pass(QualityEvaluation::PriorityQueue* _A){
    //flips,edge removal, multi-face removal operations
    bool changed = false;
    // std::cout << "Topological pass" << std::endl;
    return changed;
}
void MainLoop::edge_contraction_pass(QualityEvaluation::PriorityQueue* _A){
    // E <- edges of the triangles in A
    // foreach e in E
    // if e still exists attempt to contract e and smooth vertex created
    // return surviving set in A and the triangles in M altered by the contractions
    // std::cout << "Edge contraction pass" << std::endl;

}

void MainLoop::insertion_pass(QualityEvaluation::PriorityQueue* _A){
    // std::cout << "Insertion pass" << std::endl;

}

void MainLoop::smoothing_pass(QualityEvaluation::PriorityQueue* _A){
    // V <- the vertices in A
    std::vector<OpenMesh::SmartVertexHandle> V;
    for (unsigned long i = 0; i < _A->size(); i++) {
        QualityEvaluation::Triangle top = _A->top();
        _A->pop();
        auto fh = mesh_.face_handle(top.face_id_);
        for(auto fv_it = mesh_.fv_iter(fh); fv_it.is_valid(); ++fv_it){
            V.emplace_back(*fv_it);
        }
    }

    for(auto v: V){
        smoother_.smooth(v);
        // A <- A U the triangles adjoining v in M
        for(auto vf_it = mesh_.vf_cwiter(v); vf_it.is_valid(); ++vf_it){
            double quality = qe_.evaluate(*vf_it);
            _A->push(QualityEvaluation::Triangle(v.idx(),quality));
        }
    }
    // foreach v in V
    // Smooth v
    //std::cout << "Smoothing pass" << std::endl;
}
