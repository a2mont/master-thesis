#include "TetLoop.hh"

void TetLoop::loop(int _max_iter){
    for(int i=0; i < _max_iter; ++i){
        PriorityQueue worstTets;

        computeQuality();

        Tet currentWorst = quality_queue_.top();
        while(currentWorst.quality_ < q_min_ && !quality_queue_.empty()){
            worstTets.push(currentWorst);
            quality_queue_.pop();
            currentWorst = quality_queue_.top();
        }

        if(worstTets.empty()){
            std::cout << "\033[1;32mAll tets are of good enough quality\nWorst tet: \033[0m" << currentWorst.toString() << std::endl;
        }
        else{
            std::cout << worstTets.size() << " tets of bad quality"
                      << "\nWorst tet: " << worstTets.top().toString()
                      << std::endl;
            if(includeLogs_)
                log(logger_);
            improve_mesh(worstTets);
        }
        if(includeLogs_)
            log(logger_, true);
        reset_queue(quality_queue_);
    }
    if(includeLogs_)
        log(timeStepLogger_, true);

}

void TetLoop::computeQuality(){
    reset_queue(quality_queue_);
    for(auto c_it = mesh_.cells_begin(); c_it != mesh_.cells_end(); ++c_it){
        double quality = QualityEvaluation::evaluate(*c_it, mesh_);
        quality_queue_.push(Tet(*c_it, quality));
    }
}
void TetLoop::log(Logger* _logger, bool _endOfLine){
    computeQuality();
    _logger->logQuality(quality_queue_.top().quality_);
    if(_endOfLine){
        logger_->nextLine();
    }
}

void TetLoop::reset_queue(TetLoop::PriorityQueue& _queue){
    PriorityQueue empty;
    std::swap(_queue, empty);
}

void TetLoop::improve_mesh(PriorityQueue _badTets){
    // B <- set of triangles in M with q < q_min
    // foreach t in B
    // if t still exists and q < q_min improve_tet(M,t,q_min)
    while (!_badTets.empty()) {
        auto t = _badTets.top();
        _badTets.pop();
        if(!mesh_.is_deleted(t.cell_handle_)){
            improve_tet(t);
        }
    }
}

void TetLoop::improve_tet(Tet _t){
    // A <- t, a triangle from the mesh
    PriorityQueue A;
    A.push(_t);
    bool changed;
    for (int i = 0; i < 1; ++i) {
        int maxIter = 10;
        changed = false;
        do {
            // A <- topological_pass(A,M)
//            changed = topologial_pass(&A);
            if(A.top().quality_ >= q_min_){ //|| A.empty()){
                if(includeLogs_)
                    log(logger_,true);
                return;
            }
        } while (changed && maxIter-- > 0);

        if(includeLogs_)
            log(logger_);
        // A <- edge_contraction_pass(A,M)
        edge_contraction_pass(&A);
        if(A.top().quality_ >= q_min_ || A.empty()){
            if(includeLogs_)
                log(logger_, true);
            return;
        }
        if(includeLogs_)
            log(logger_);
        // A <- insertion_pass(A,M)
//        insertion_pass(&A);
        if(A.top().quality_ >= q_min_){ //|| A.empty()){
            if(includeLogs_)
                log(logger_, true);
            return;
        }
        if(includeLogs_)
            log(logger_);
        // smoothing_pass()
        smoothing_pass(&A, 1);
        if(A.top().quality_ >= q_min_ || A.empty()){
            if(includeLogs_)
                log(logger_, true);
            return;
        }
        if(includeLogs_)
            log(logger_);
    }

}

bool TetLoop::topologial_pass(PriorityQueue* _A){
    //flips,edge removal, multi-face removal operations
    std::cout << "\033[1;36mTopological pass with _A of size: \033[0m" <<  _A->size() << std::endl;
    bool changed = false;
    // for t in A that still exists
    //      for each edge e of t (if t still exists)
    //          attempt to remove e
    //      for each f of t (if t still exists)
    //          attempt to remove f (multi-face, 2-3/2-2 flip)
    // A <- A U new tets created by the operation
    while(!_A->empty()){
        auto top = _A->top();
        _A->pop();
        auto ch = top.cell_handle_;
        std::vector<EdgeHandle> edges;
        std::vector<FaceHandle> faces;
        for(auto c_eh: mesh_.cell_edges(ch)){
            edges.push_back(c_eh);
        }
        for(auto c_fh: mesh_.cell_faces(ch)){
            faces.push_back(c_fh);
        }
        for(auto e: edges){
            if(!mesh_.is_deleted(e)){
                // remove e ?
            }
        }
        for(auto f: faces){
            if(!mesh_.is_deleted(f)){
                // remove f ?
            }
        }

    }

    return changed;
}
void TetLoop::edge_contraction_pass(PriorityQueue* _A){
    std::cout << "\033[1;36mEdge contraction pass with _A of size: \033[0m" <<  _A->size() << std::endl;
    std::vector<EdgeHandle> E;
    std::vector<CellHandle> newTets;
    while(!_A->empty()){
        E.clear();
        auto top = _A->top();
        _A->pop();
        auto ch = top.cell_handle_;
        // E <- edges of the tet in A
        for(auto eh: mesh_.cell_edges(ch)){
            E.push_back(eh);
        }
        for(auto e: E){
            if(!mesh_.is_deleted(e)){
                auto he = mesh_.edge_halfedges(e)[0];
                auto from = mesh_.from_vertex_handle(he);
                auto to = mesh_.to_vertex_handle(he);
                std::cout << "Boundary: "<< mesh_.is_boundary(from) << " " << mesh_.is_boundary(to) << std::endl;
//                if(mesh_.is_boundary(from) && mesh_.is_boundary(to))
//                    continue;
                // if only one vertex is on the boundary, invert the collapse direction
                if(mesh_.is_boundary(from)){
                    he = mesh_.edge_halfedges(e)[1];
                    from = mesh_.from_vertex_handle(he);
                    to = mesh_.to_vertex_handle(he);
                }
                // get the tets altered by contraction
                for(auto v_ch: mesh_.vertex_cells(from)){
                    newTets.emplace_back(v_ch);
                }
                std::cout << "Constraint vertices: " << from.idx() << " "  << to.idx() <<" collapse" << std::endl;
                if(constraint_vhs_.count(from.idx()) > 0){
                    constraint_vhs_.erase(from.idx());
                    constraint_vhs_[to.idx()] = to.idx();
                }
                auto remain = mesh_.collapse_edge(he);
                Smoothing::smooth(mesh_, remain);
//                break;
            }
        }

    }
    // return surviving set in A and the tets in M altered by the contractions
    for(auto new_fh: newTets){
        if (!mesh_.is_deleted(new_fh)) {
            _A->push(Tet(new_fh, QualityEvaluation::evaluate(new_fh,mesh_)));
        }
    }

}



void TetLoop::insertion_pass(PriorityQueue* _A){
    std::cout << "\033[1;36mInsertion pass with _A of size: \033[0m" <<  _A->size()  << std::endl;


}

void TetLoop::smoothing_pass(PriorityQueue* _A, int _iterations){
    std::cout << "\033[1;36mSmoothing pass with _A of size: \033[0m" <<  _A->size() << std::endl;
    // V the vertices of the tets in A
    for(int i = 0; i < _iterations; ++i){
        std::vector<OpenVolumeMesh::VertexHandle> V;
        while(!_A->empty()) {
            Tet top = _A->top();
            _A->pop();
            auto ch = top.cell_handle_;
            for(auto cvh: mesh_.cell_vertices(ch)){
                V.emplace_back(cvh);
            }
        }
        for(auto v: V){
            Smoothing::smooth(mesh_, v);
            for(auto vch: mesh_.vertex_cells(v)){
                double quality = QualityEvaluation::evaluate(vch, mesh_);
                _A->push(Tet(vch,quality));
            }
        }
    }
}
