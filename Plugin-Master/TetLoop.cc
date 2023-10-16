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
void TetLoop::computeQuality(TetLoop::PriorityQueue& _queue, TetrahedralMesh& _mesh){
    reset_queue(_queue);
    for(auto c_it = _mesh.cells_begin(); c_it != _mesh.cells_end(); ++c_it){
        double quality = QualityEvaluation::evaluate(*c_it, _mesh);
        _queue.push(Tet(*c_it, quality));
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
//        edge_contraction_pass(&A);
        if(A.top().quality_ >= q_min_) {//|| A.empty()){
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
//        smoothing_pass(&A, 1);
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
    bool changed, edge, face = false;
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
        if(mesh_.is_deleted(ch))
            continue;

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
                edge = edgeRemoval(e);
                changed = changed || edge;
            }
        }
        for(auto f: faces){
            if(!mesh_.is_deleted(f)){
                face = faceRemoval(f);
                changed = changed || face;
            }
        }
    }

    return changed;
}

bool TetLoop::edgeRemoval(EdgeHandle _eh, bool _verbose){
    bool changed = false;
    // let I the set of tets including eh
    std::set<CellHandle> tetsIncludingEdge;
    // let R the edge ring around eh
    std::set<VertexHandle> oneRingVertices;
    std::vector<EdgeHandle> newEdges;
    PriorityQueue localQueue;

    auto from = mesh_.edge(_eh).from_vertex();
    auto to = mesh_.edge(_eh).to_vertex();

    double q_old = 0;
    double maxQuality = 0;
    int maxId = -1;
    int current = 0;

    for(auto cell: mesh_.edge_cells(_eh)){
        tetsIncludingEdge.insert(cell);
        for(auto vh: mesh_.cell_vertices(cell)){
            if(vh.idx() != from.idx() && vh.idx() != to.idx()){
                oneRingVertices.insert(vh);
            }
        }
    }

    // Empirical results from paper, edge removal for m > 7 rarely improves quality
    if(oneRingVertices.size() > 7)
        return changed = false;

    for(auto ch: tetsIncludingEdge){
        double quality = QualityEvaluation::evaluate(ch, mesh_);
        localQueue.push(Tet(ch,quality));
    }

    q_old = localQueue.top().quality_;

    reset_queue(localQueue);

    auto tempBaseMesh = mesh_;
    auto newVertex = tempBaseMesh.split_edge(_eh);

    for(auto vheh: tempBaseMesh.outgoing_halfedges(newVertex)){
        if(tempBaseMesh.to_vertex_handle(vheh) == to || tempBaseMesh.to_vertex_handle(vheh) == from)
            continue;
        newEdges.push_back(tempBaseMesh.edge_handle(vheh));
    }

    for(auto eh: newEdges){
        auto tempMesh = tempBaseMesh;
        auto hehs = tempMesh.edge_halfedges(eh);
        auto heh = tempMesh.from_vertex_handle(hehs[0]) == newVertex ? hehs[0] : hehs[1];
        if(!link_condition(tempMesh, heh))
            continue;
        auto tempNew = tempMesh.collapse_edge(heh);
        Smoothing::smooth(tempMesh, tempNew);
        for(auto ch: tempMesh.vertex_cells(tempNew)){
            if(!tempMesh.is_deleted(ch)){
                double quality = QualityEvaluation::evaluate(ch, tempMesh);
                localQueue.push(Tet(ch,quality));
            }
        }
        if(!localQueue.empty() && localQueue.top().quality_ > maxQuality){
            maxQuality = localQueue.top().quality_;
            maxId = current;
        }
        reset_queue(localQueue);
        current++;
    }

    if(maxQuality > q_old){
        auto edgeToCollapse = newEdges[maxId];
        auto hehs = tempBaseMesh.edge_halfedges(edgeToCollapse);
        auto heh = tempBaseMesh.from_vertex_handle(hehs[0]) == newVertex ? hehs[0] : hehs[1];
        if(link_condition(tempBaseMesh, heh)){
            auto newFrom = tempBaseMesh.collapse_edge(heh);
            mesh_ = tempBaseMesh;
            Smoothing::smooth(mesh_, newFrom);
            return changed = true;
        }
    }
    if(_verbose){
        std::cout
                << "\033[1;33mEdge removal reverted, old: \033[0m"
                << q_old
                << "\033[1;33m vs max: \033[0m"
                << maxQuality
                << std::endl;
    }

    return changed;
}

bool TetLoop::faceRemoval(FaceHandle _fh){
    bool changed = false;
    auto tempMesh = mesh_;
    auto bestMesh = mesh_;
    PriorityQueue tempQueue;

    computeQuality(tempQueue, mesh_);
    double q_old = tempQueue.top().quality_;
    double q_max = q_old;

    flip23(_fh);
    double q_new = tempQueue.top().quality_;
    computeQuality(tempQueue, mesh_);
    q_new = tempQueue.top().quality_;
    if(q_new > q_max){
        q_max = q_new;
        bestMesh = mesh_;
        changed = true;
//        std::cout << "2-3 flip has better quality" << std::endl;
    }

    mesh_ = tempMesh;
    flip22(_fh);
    computeQuality(tempQueue, mesh_);
    if(q_new > q_max){
        q_max = q_new;
        bestMesh = mesh_;
        changed = true;
//        std::cout << "2-2 flip has better quality" << std::endl;
    }

//    mesh_ = tempMesh;
//    multiFace(_fh);
//    computeQuality(tempQueue, mesh_);
//    if(q_new > q_max){
//        q_max = q_new;
//        bestMesh = mesh_;
//        changed = true;
//        std::cout << "Multi face has better quality" << std::endl;
//    }

    mesh_ = bestMesh;

    return changed;
}

void TetLoop::multiFace(FaceHandle _fh){
    if(mesh_.is_boundary(_fh))
        return;
    auto hfhs = mesh_.face_halffaces(_fh);
    auto a = mesh_.halfface_opposite_vertex(hfhs[0]);
    auto b = mesh_.halfface_opposite_vertex(hfhs[1]);

    std::vector<VertexHandle> faceVertices; // u,v,w

    for(auto vh: mesh_.get_halfface_vertices(hfhs[0])){
        faceVertices.push_back(vh);
    }
    auto results_uv = testNeighbor(a,b,faceVertices[0],faceVertices[1]);
    auto results_vw = testNeighbor(a,b,faceVertices[1],faceVertices[2]);
    auto results_wu = testNeighbor(a,b,faceVertices[2],faceVertices[0]);

    double q_old = std::min({QualityEvaluation::evaluate(a,faceVertices[0],faceVertices[1],faceVertices[2], mesh_),
                             QualityEvaluation::evaluate(faceVertices[0],faceVertices[1],faceVertices[2],b, mesh_),
                             results_uv.o_,
                             results_vw.o_,
                             results_wu.o_});
    double q_new = std::min({results_uv.n_, results_vw.n_, results_wu.n_});
    if(q_new > q_old){
        flip23(_fh);
        std::vector<FaceWithChildren> total;
        total.insert(total.end(), results_uv.h_.begin(), results_uv.h_.end());
        total.insert(total.end(), results_vw.h_.begin(), results_vw.h_.end());
        total.insert(total.end(), results_wu.h_.begin(), results_wu.h_.end());

        for(auto g: total){
            flip32Recurse(g, FaceWithChildren(_fh));
        }
    }
}

/**
 * With f the face to remove
 * @param a, a vertex sandwiching f
 * @param b, a vertex sandwiching f
 * @param u, vertex of f
 * @param w, vertex of f
 */
TetLoop::TestNeighborResult TetLoop::testNeighbor(VertexHandle a,
                           VertexHandle b,
                           VertexHandle u,
                           VertexHandle w){
    double q_uw = QualityEvaluation::evaluate(a,b,u,w,mesh_);
    auto he_uw = mesh_.halfedge(u,w);
    auto uw = mesh_.edge_handle(he_uw);
    FaceWithChildren g;
    VertexHandle v;
    bool correctBoundary = !mesh_.is_boundary(uw) ||
            (mesh_.is_boundary(a) && mesh_.is_boundary(b));

    if(mesh_.valence(uw) == 4 && correctBoundary){
        for(auto hehfh: mesh_.halfedge_halffaces(he_uw)){
            auto next_he = mesh_.next_halfedge_in_halfface(he_uw, hehfh);
            if(mesh_.to_vertex_handle(next_he) == a || mesh_.to_vertex_handle(next_he) == b){
                continue;
            }
            g = mesh_.face_handle(hehfh);
            v = mesh_.to_vertex_handle(next_he);
            break;
        }
        double j_uv = orient3D(a,b,u,v);
        double j_vw = orient3D(a,b,v,w);
        double j_wu = orient3D(a,b,w,u);
        if((j_uv > 0 && j_vw > 0) ||
           (j_vw > 0 && j_wu > 0) ||
           (j_wu > 0 && j_uv > 0)){
            auto result_uv = testNeighbor(a,b,u,v);
            auto result_vw = testNeighbor(a,b,v,w);
            if(!result_uv.h_.empty()){
                for(auto h: result_uv.h_){
                    g.children_.insert(h.fh_);
                }
            }
            if(!result_vw.h_.empty()){
                for(auto h: result_vw.h_){
                    g.children_.insert(h.fh_);
                }
            }
            double q_old = std::min({QualityEvaluation::evaluate(a,u,v,w,mesh_),
                                     QualityEvaluation::evaluate(u,v,w,b,mesh_),
                                     result_uv.o_,
                                     result_vw.o_});
            double q_new = std::min({result_uv.n_, result_vw.n_});
            if(q_new > q_old || q_new > q_uw){
                return TestNeighborResult(q_old, q_new, std::vector<FaceWithChildren> {g});
            }
        }
    }

    return TestNeighborResult(QualityEvaluation::getMaxQuality(), q_uw, std::vector<FaceWithChildren> {});
}

double TetLoop::orient3D(VertexHandle _a, VertexHandle _b, VertexHandle _c, VertexHandle _d){
    std::vector<Point> points = {mesh_.vertex(_a), mesh_.vertex(_b), mesh_.vertex(_c), mesh_.vertex(_d)};
    std::vector<Eigen::Vector3d> cols;
    for(int i = 0; i < 3; ++i){
        cols.push_back(Eigen::Vector3d(
                           points[0][i] - points[3][i],
                           points[1][i] - points[3][i],
                           points[2][i] - points[3][i]));
    }
    Eigen::Matrix3d matrix;
    int j = 0;
    for(auto c: cols){
        matrix.col(j++) = c;
    }
    return matrix.determinant();
}

void TetLoop::flip32Recurse(TetLoop::FaceWithChildren _g, TetLoop::FaceWithChildren _parent){
    EdgeHandle toRemove;
    for(auto parentEdge: mesh_.face_edges(_parent.fh_)){
        for(auto edge: mesh_.face_edges(_g.fh_)){
            if(edge == parentEdge){
                toRemove = edge;
            }
        }
    }
    flip32(toRemove);
    for(auto h: _g.children_){
        flip32Recurse(h, _g);
    }
}

void TetLoop::flip32(EdgeHandle _eh){
    if(mesh_.is_deleted(_eh))
        return;
    auto ends = mesh_.edge_vertices(_eh);
    HalfEdgeHandle toCollapse = mesh_.InvalidHalfEdgeHandle;
    auto temp = mesh_;
    auto newVertex = temp.split_edge(_eh);
    for(auto heh: mesh_.outgoing_halfedges(newVertex)){
        if(mesh_.to_vertex_handle(heh) != ends[0] && mesh_.to_vertex_handle(heh) != ends[1]){
            toCollapse = heh;
        }
    }
    if(link_condition(temp,toCollapse)){
        temp.collapse_edge(toCollapse);
        mesh_ = temp;
    }
}

void TetLoop::flip23(FaceHandle _fh){
    if(mesh_.is_deleted(_fh))
        return;
    auto hfh = mesh_.face_halffaces(_fh)[0];
    auto toVertex = mesh_.halfface_opposite_vertex(hfh);
    if(mesh_.is_boundary(toVertex)){
        hfh = mesh_.face_halffaces(_fh)[1];
        toVertex = mesh_.halfface_opposite_vertex(hfh);
    }
    auto temp = mesh_;
    auto newVertex = temp.split_face(_fh);
    auto heToCollapse = temp.halfedge(newVertex, toVertex);
    if(link_condition(temp, heToCollapse)){
        temp.collapse_edge(heToCollapse);
        mesh_ = temp;
    }
}

void TetLoop::flip22(FaceHandle _fh){
    if(mesh_.is_deleted(_fh))
        return;
    HalfEdgeHandle heh;
    auto hfh = mesh_.face_halffaces(_fh)[0];
    for(auto he: mesh_.face_halfedges(_fh)){
        heh = he;
    }
    auto temp = mesh_;
    auto toVertex = temp.to_vertex_handle(mesh_.next_halfedge_in_halfface(heh, hfh));
    auto newVertex = temp.split_edge(temp.edge_handle(heh));
    auto heToCollapse = temp.halfedge(newVertex, toVertex);
    if(link_condition(temp, heToCollapse)){
        temp.collapse_edge(heToCollapse);
        mesh_ = temp;
    }
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
            auto remain = contractEdge(e, newTets);
            if(remain.is_valid())
                Smoothing::smooth(mesh_, remain);
        }

    }
    // return surviving set in A and the tets in M altered by the contractions
    for(auto new_fh: newTets){
        if (!mesh_.is_deleted(new_fh)) {
            _A->push(Tet(new_fh, QualityEvaluation::evaluate(new_fh,mesh_)));
        }
    }

}

VertexHandle TetLoop::contractEdge(EdgeHandle _eh, std::vector<CellHandle>& _tetsAltered){
    if(mesh_.is_deleted(_eh)){
        return mesh_.InvalidVertexHandle;
    }
    auto he = mesh_.edge_halfedges(_eh)[0];
    auto from = mesh_.from_vertex_handle(he);
    auto to = mesh_.to_vertex_handle(he);
    if((mesh_.is_boundary(from) && mesh_.is_boundary(to))|| !link_condition(mesh_, _eh))
        return mesh_.InvalidVertexHandle;
    // if only one vertex is on the boundary, invert the collapse direction
    if(mesh_.is_boundary(from)){
        he = mesh_.edge_halfedges(_eh)[1];
        from = mesh_.from_vertex_handle(he);
        to = mesh_.to_vertex_handle(he);
    }
    // get the tets altered by contraction
    for(auto v_ch: mesh_.vertex_cells(from)){
        _tetsAltered.emplace_back(v_ch);
    }
    std::cout << "Constraint vertices: " << from.idx() << " "  << to.idx() <<" collapse" << std::endl;
    if(constraint_vhs_.count(from.idx()) > 0){
        constraint_vhs_.erase(from.idx());
        constraint_vhs_[to.idx()] = to.idx();
    }
    auto remain = mesh_.collapse_edge(he);
    return remain;
}


void TetLoop::insertion_pass(PriorityQueue* _A){
    std::cout << "\033[1;36mInsertion pass with _A of size: \033[0m" <<  _A->size()  << std::endl;


}

void TetLoop::smoothing_pass(PriorityQueue* _A, int _iterations){
    std::cout << "\033[1;36mSmoothing pass with _A of size: \033[0m" <<  _A->size() << std::endl;
    // V the vertices of the tets in A
    for(int i = 0; i < _iterations; ++i){
        std::vector<VertexHandle> V;
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
