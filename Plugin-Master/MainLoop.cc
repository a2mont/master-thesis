#include "MainLoop.hh"



void MainLoop::loop(ACG::Vec3d _displacement, int _constraint_vh, bool _verbose,int _max_iter ){

    for(int i=0; i < _max_iter; ++i){
        PriorityQueue worstTriangles;
        Triangle worstTriangle(0,0);

        vd_.displace(_displacement, _constraint_vh, _verbose);

        for(auto fh: mesh_.faces()){
            double quality = QualityEvaluation::evaluate(fh, mesh_,_verbose);
            quality_queue_.push(Triangle(fh.idx(), quality));
        }

        worstTriangle = quality_queue_.top();
        while(worstTriangle.quality_ < q_min_){
            worstTriangles.push(worstTriangle);
            quality_queue_.pop();
            worstTriangle = quality_queue_.top();
        }
        if(worstTriangles.empty())
            std::cout << "All triangles are of good enough quality\nWorst triangle: " << worstTriangle.toString() << std::endl;
        else
            std::cout << worstTriangles.size() << " triangles of bad quality"
                      << "\nWorst triangle: " << worstTriangle.toString()
                      << std::endl;

        improve_mesh(worstTriangles);

        reset_queue(quality_queue_);
        mesh_.garbage_collection();
    }
}

void MainLoop::reset_queue(MainLoop::PriorityQueue& _queue){
    PriorityQueue empty;
    std::swap(_queue, empty);
}

void MainLoop::improve_mesh(PriorityQueue _badTriangles){
    // B <- set of triangles in M with q < q_min
    // foreach t in B
    // if t still exists and q < q_min improve_tet(M,t,q_min)
    while (!_badTriangles.empty()) {
        auto t = _badTriangles.top();
        _badTriangles.pop();
        // TODO check if t exists
        if(!mesh_.status(mesh_.face_handle(t.face_id_)).deleted())
            improve_triangle(t);
    }
}

void MainLoop::improve_triangle(Triangle _t){
    // A <- t, a triangle from the mesh
    PriorityQueue A;
    A.push(_t);
    for (int i = 0; i < 10; ++i) {
        bool changed = false;
        do {
            // A <- topological_pass(A,M)
            //changed = topologial_pass(&A);
            if(A.top().quality_ >= q_min_)
                return;
        } while (changed);

        // A <- edge_contraction_pass(A,M)
        edge_contraction_pass(&A);
        if(A.top().quality_ >= q_min_)
            return;

        // A <- insertion_pass(A,M)
//        insertion_pass(&A);
        if(A.top().quality_ >= q_min_)
            return;

        // smoothing_pass()
        smoothing_pass(&A);
        if(A.top().quality_ >= q_min_)
            return;
    }
}

bool MainLoop::topologial_pass(PriorityQueue* _A){
    //flips,edge removal, multi-face removal operations
    std::cout << "Topological pass" << std::endl;
    bool changed = false;
    // for each tri in a
    // for each edge in tri
    //  attempt to remove e
    for(unsigned long i=0; i < _A->size(); ++i){
        auto top = _A->top();
        _A->pop();
        auto fh = mesh_.face_handle(top.face_id_);
        for (auto fe_it = mesh_.fe_cwiter(fh); fe_it.is_valid(); ++fe_it) {
            if(!mesh_.status(fh).deleted()){
                mesh_.delete_edge(*fe_it);
                changed = true;
            }
        }
    }

    return changed;
}
void MainLoop::edge_contraction_pass(PriorityQueue* _A){
//    std::cout << "Edge contraction pass" << std::endl;
    // E <- edges of the triangles in A
    std::vector<OpenMesh::SmartEdgeHandle> E;
    std::vector<Triangle> newTriangles;
    for(unsigned long i=0; i < _A->size(); ++i){
        auto top = _A->top();
        _A->pop();
        auto fh = mesh_.face_handle(top.face_id_);

        for (auto fe_it = mesh_.fe_cwiter(fh); fe_it.is_valid(); ++fe_it) {
                E.emplace_back(*fe_it);
        }
        // foreach e in E
        // if e still exists attempt to contract e and smooth vertex created
        for(auto e: E){
            if(!mesh_.status(fh).deleted()){
                if (mesh_.is_collapse_ok(e.h0())) {
                    auto he = e.h0();
                    // get the faces altered by contraction
                    std::vector<OpenMesh::SmartFaceHandle> nextFaces {e.h0().face(), e.h1().face()};
                    for(auto fh: nextFaces){
                        newTriangles.emplace_back(Triangle(fh.idx(), QualityEvaluation::evaluate(fh, mesh_)));
                    }
                    mesh_.collapse(he);
                    smoother_.smooth(he.to());
                    break;
                }
            }
        }
        E.clear();
    }

    // return surviving set in A and the triangles in M altered by the contractions
    for(auto triangle: newTriangles){
        _A->push(triangle);
    }

}



void MainLoop::insertion_pass(PriorityQueue* _A){
    // std::cout << "Insertion pass" << std::endl;

}

void MainLoop::smoothing_pass(PriorityQueue* _A){
//    std::cout << "Smoothing pass" << std::endl;
    // V <- the vertices in A
    std::vector<OpenMesh::SmartVertexHandle> V;
    for (unsigned long i = 0; i < _A->size(); ++i) {
        Triangle top = _A->top();
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
            double quality = QualityEvaluation::evaluate(*vf_it, mesh_);
            _A->push(Triangle(vf_it->idx(),quality));
        }
    }
}
