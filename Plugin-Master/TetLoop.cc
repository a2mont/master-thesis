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
//    computeQuality();
//    while (!quality_queue_.empty()) {
//        auto tet = quality_queue_.top();
//        quality_queue_.pop();
//        if(tet.quality_ < 10e-6){
//            mesh_.delete_cell(tet.cell_handle_);
//        }
//    }
    if(includeLogs_)
        log(timeStepLogger_, true);

}

void TetLoop::computeQuality(){
    reset_queue(quality_queue_);
    for(auto c_it = mesh_.cells_begin(); c_it != mesh_.cells_end(); ++c_it){
        if(!c_it->is_valid() || mesh_.is_deleted(*c_it)){
            std::cout <<*c_it <<" is deleted" << std::endl;
            continue;
        }
        double quality = QualityEvaluation::evaluate(*c_it, mesh_);
//        std::cout << quality << std::endl;
        quality_queue_.push(Tet(*c_it, quality));
//        std::cout << "after" << std::endl;
    }
//    std::cout << "end" << std::endl;
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
    std::cout << "Improve tet: "<< _t.cell_handle_ << std::endl;
    for (int i = 0; i < 1; ++i) {
        int maxIter = 10;
        changed = false;
        std::cout << "\033[1;36mTopological pass\033[0m"<< std::endl;
        do {
            // A <- topological_pass(A,M)
//            changed = topologial_pass(&A);
            if(A.top().quality_ >= q_min_ || A.empty()){
                if(includeLogs_)
                    log(logger_,true);
                return;
            }
        } while (changed && maxIter-- > 0);

        if(includeLogs_)
            log(logger_);
        // A <- edge_contraction_pass(A,M)
        std::cout << "\033[1;36mEdge contraction pass \033[0m" << std::endl;
//        edge_contraction_pass(&A);
        if(A.top().quality_ >= q_min_) {//|| A.empty()){
            if(includeLogs_)
                log(logger_, true);
            return;
        }
        if(includeLogs_)
            log(logger_);
        // A <- insertion_pass(A,M)
        std::cout << "\033[1;36mInsertion pass\033[0m" << std::endl;
        insertion_pass(&A);
        if(A.top().quality_ >= q_min_){ //|| A.empty()){
            if(includeLogs_)
                log(logger_, true);
            return;
        }
        if(includeLogs_)
            log(logger_);
        // smoothing_pass()
        std::cout << "\033[1;36mSmoothing pass\033[0m" << std::endl;
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

// ----------- ** Topological pass ** ---------------------

bool TetLoop::topologial_pass(PriorityQueue* _A){
    //flips,edge removal, multi-face removal operations
    bool changed, edge, face;
    changed = edge = face = false;
    std::vector<CellHandle> cellsAdded;
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
        for(auto e: edges){
            if(!mesh_.is_deleted(e)){
                edge = edgeRemoval(e, &cellsAdded,false);
                changed = changed || edge;
            }
        }
        if(mesh_.is_deleted(ch))
            continue;
        for(auto c_fh: mesh_.cell_faces(ch)){
            faces.push_back(c_fh);
        }
        for(auto f: faces){
            if(!mesh_.is_deleted(f)){
                face = faceRemoval(f, false);
                changed = changed || face;
            }
        }
    }
//    std::cout << "Tets added: "<< cellsAdded.size() << std::endl;

    return changed;
}

// Call when we don't track the cells added by the process (e.g. Tests.cc)
bool TetLoop::edgeRemoval(EdgeHandle _eh, bool _verbose){
    std::vector<CellHandle> ignore;
    return edgeRemoval(_eh, &ignore, _verbose);
}

bool TetLoop::edgeRemoval(EdgeHandle _eh, std::vector<CellHandle>* const _cellsAdded, bool _verbose){
    bool changed = false;
    if(mesh_.is_deleted(_eh) || mesh_.is_boundary(_eh))
        return changed;

    // let I the set of tets including eh
    std::set<CellHandle> tetsIncludingEdge;
    // let R the edge ring around eh
    std::set<VertexHandle> oneRingVertices;
    std::vector<HalfEdgeHandle> newHalfEdges;
    PriorityQueue localQueue;

    HalfEdgeHandle toCollapse(-1);

    auto from = mesh_.edge(_eh).from_vertex();
    auto to = mesh_.edge(_eh).to_vertex();

    double q_old = 0;
    double maxQuality = 0;

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
    VertexHandle newVertex = tempBaseMesh.split_edge(_eh);

    for(auto vheh: tempBaseMesh.outgoing_halfedges(newVertex)){
        if(tempBaseMesh.to_vertex_handle(vheh) == to || tempBaseMesh.to_vertex_handle(vheh) == from)
            continue;
        newHalfEdges.push_back(vheh);
    }
    for(auto heh: newHalfEdges){
        auto tempMesh = tempBaseMesh;
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
            toCollapse = heh;
        }
        reset_queue(localQueue);
    }

    if(maxQuality > q_old){
        if(link_condition(tempBaseMesh, toCollapse) && toCollapse.is_valid()){
            auto addedVertex = mesh_.split_edge(_eh);
            std::vector<CellHandle> adjacentCells;
            for(auto ch: mesh_.halfedge_cells(toCollapse)){
                adjacentCells.push_back(ch);
            }
            for(auto ch: mesh_.vertex_cells(addedVertex)){
                if(std::find(adjacentCells.begin(), adjacentCells.end(), ch) == adjacentCells.end()){
                    _cellsAdded->push_back(ch);
                }
            }
            auto newFrom = mesh_.collapse_edge(toCollapse);
            Smoothing::smooth(mesh_, newFrom);
            changed = true;
            return changed;
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

bool TetLoop::faceRemoval(FaceHandle _fh, bool _verbose){
    bool changed = false;
    auto tempMesh = mesh_;
    TetrahedralMesh bestMesh = mesh_;
    PriorityQueue tempQueue;

    /* track the changes of the best remeshing
     * ids:
     *  - 0: 2-3 split
     *  - 1: 2-2 split
     *  - 2: multiface
    */
    std::vector<std::vector<CellHandle>> addedTets(3);
    int bestId = -1;

    computeQuality(tempQueue, mesh_);
    double q_old = tempQueue.top().quality_;
    double q_max = q_old;
    double q_new = 0;

//    flip23(_fh, &addedTets[0]);
//    computeQuality(tempQueue, mesh_);
//    q_new = tempQueue.top().quality_;
////    std::cout << "2-3\nNew: " << q_new << " Max: "<< q_max << std::endl;
//    if(q_new > q_max){
//        q_max = q_new;
//        bestMesh = mesh_;
//        changed = true;
//        bestId = 0;
////        std::cout << "2-3 flip has better quality" << std::endl;
////        std::cout << "2-3 added tets: "<< addedTets[0].size() << std::endl;
//    }

//    mesh_ = tempMesh;
//    flip22(_fh, &addedTets[1]);
//    computeQuality(tempQueue, mesh_);
//    q_new = tempQueue.top().quality_;
////    std::cout << "2-2\nNew: " << q_new << " Max: "<< q_max << std::endl;
//    if(q_new > q_max){
//        q_max = q_new;
//        bestMesh = mesh_;
//        changed = true;
//        bestId = 1;
////        std::cout << "2-2 flip has better quality" << std::endl;
////        std::cout << "2-2 added tets: "<< addedTets[1].size() << std::endl;
//    }

//    mesh_ = tempMesh;
    multiFace(_fh, &addedTets[2]);
    computeQuality(tempQueue, mesh_);
    q_new = tempQueue.top().quality_;
//    std::cout << "Multi\nNew: " << q_new << " Max: "<< q_max << std::endl;
    if(q_new > q_max){
        q_max = q_new;
        bestMesh = mesh_;
        changed = true;
        bestId = 2;
//        std::cout << "Multi face has better quality" << std::endl;
//        std::cout << "multi added tets: "<< addedTets[2].size() << std::endl;
    }

    mesh_ = bestMesh;
    std::cout << "---------------------------------------" << std::endl;
    if(bestId != -1){
        std::cout << "Best remeshing with id: " << bestId << ", #added tets: "<< addedTets[bestId].size() << std::endl;
    }

    if(_verbose && !changed){
        std::cout
                << "\033[1;33mFace removal reverted, old: \033[0m"
                << q_old
                << "\033[1;33m vs max: \033[0m"
                << q_max
                << std::endl;
    }

    return changed;
}

// Implementation from paper Two Discrete Optimization Algorithms for the Topological Improvement of
// Tetrahedral Meshes
void TetLoop::multiFace(FaceHandle _fh, std::vector<CellHandle>* const _cellsAdded){
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
                             results_wu.o_
                            });
    double q_new = std::min({results_uv.n_, results_vw.n_, results_wu.n_});
    if(q_new > q_old){
        std::cout << "Dans multi: "<< q_old << ", vs q_new: "<< q_new << std::endl;
        flip23(_fh, _cellsAdded);
        std::vector<FaceWithChildren> total;
        total.insert(total.end(), results_uv.h_.begin(), results_uv.h_.end());
        total.insert(total.end(), results_vw.h_.begin(), results_vw.h_.end());
        total.insert(total.end(), results_wu.h_.begin(), results_wu.h_.end());
        for(auto g: total){
            flip32Recurse(g, FaceWithChildren(_fh), _cellsAdded);
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

void TetLoop::flip32Recurse(TetLoop::FaceWithChildren _g,
                            TetLoop::FaceWithChildren _parent,
                            std::vector<CellHandle>* const _cellsAdded){
    EdgeHandle toRemove;
    for(auto parentEdge: mesh_.face_edges(_parent.fh_)){
        for(auto edge: mesh_.face_edges(_g.fh_)){
            if(edge == parentEdge){
                toRemove = edge;
            }
        }
    }
    flip32(toRemove, _cellsAdded);
    for(auto h: _g.children_){
        flip32Recurse(h, _g, _cellsAdded);
    }
}

void TetLoop::flip32(EdgeHandle _eh, std::vector<CellHandle>* const _cellsAdded){
    if(mesh_.is_deleted(_eh) || mesh_.is_boundary(_eh))
        return;

    auto ends = mesh_.edge_vertices(_eh);
    HalfEdgeHandle toCollapse(-1);
    auto temp = mesh_;
    VertexHandle newVertex(-1);
    try {
       newVertex = temp.split_edge(_eh);
    } catch (...) {
        std::cout << "------------ Error splitting edge: "<< _eh << "--------------" << std::endl;
    }
    if(!newVertex.is_valid() || mesh_.is_boundary(ends[0]) || mesh_.is_boundary(ends[1])){
        return;
    }
    for(auto heh: mesh_.outgoing_halfedges(newVertex)){
        if(mesh_.to_vertex_handle(heh) != ends[0] && mesh_.to_vertex_handle(heh) != ends[1]){
            toCollapse = heh;
        }
    }
    if(link_condition(temp,toCollapse)){
        auto addedVertex = mesh_.split_edge(_eh);
        std::vector<CellHandle> adjacentCells;
        for(auto ch: mesh_.halfedge_cells(toCollapse)){
            adjacentCells.push_back(ch);
        }
        for(auto ch: mesh_.vertex_cells(addedVertex)){
            if(std::find(adjacentCells.begin(), adjacentCells.end(), ch) == adjacentCells.end()){
                _cellsAdded->push_back(ch);
            }
        }
        mesh_.collapse_edge(toCollapse);
//        std::cout << "flip 3-2" << std::endl;
    }
}

void TetLoop::flip23(FaceHandle _fh, std::vector<CellHandle>* const _cellsAdded){
    if(mesh_.is_deleted(_fh) || mesh_.is_boundary(_fh))
        return;
    auto hfh = mesh_.face_halffaces(_fh)[0];
    auto toVertex = mesh_.halfface_opposite_vertex(hfh);
    if(mesh_.is_boundary(toVertex)){
        hfh = mesh_.face_halffaces(_fh)[1];
        toVertex = mesh_.halfface_opposite_vertex(hfh);
    }
    if(mesh_.is_boundary(toVertex))
        return;

    auto hfh_opp = mesh_.opposite_halfface_handle(hfh);
    auto vh_opp = mesh_.halfface_opposite_vertex(hfh_opp);

    auto temp = mesh_;
    auto newVertex = temp.split_face(_fh);
    auto heToCollapse = temp.halfedge(newVertex, toVertex);
    auto remainingHe = temp.halfedge(vh_opp, newVertex);

    if(link_condition(temp, heToCollapse)){
        mesh_.split_face(_fh);
        for(auto ch: mesh_.halfedge_cells(remainingHe)){
            _cellsAdded->push_back(ch);
        }
        std::cout << "flip 2-3" << std::endl;
        mesh_.collapse_edge(heToCollapse);
    }
}

void TetLoop::flip22(FaceHandle _fh, std::vector<CellHandle>* const _cellsAdded){
    if(mesh_.is_deleted(_fh)|| mesh_.is_boundary(_fh))
        return;
    HalfEdgeHandle heh(-1);
    HalfEdgeHandle toCollapse(-1);
    auto hfh = mesh_.face_halffaces(_fh)[0];
    for(auto he: mesh_.face_halfedges(_fh)){
        if(!mesh_.is_boundary(he)){
            heh = he;
            break;
        }
    }
    // Get the heh going in the opposite direction from collapse (_nc = next cell)
    auto adjHf = mesh_.adjacent_halfface_in_cell(hfh, heh);
    auto oppHf_nc = mesh_.opposite_halfface_handle(adjHf);
    auto oppHe_nc = mesh_.opposite_halfedge_handle(heh);
    auto adjHf_nc = mesh_.adjacent_halfface_in_cell(oppHf_nc, oppHe_nc);
    auto oppVh_nc = mesh_.to_vertex_handle(mesh_.next_halfedge_in_halfface(oppHe_nc, adjHf_nc));

    auto temp = mesh_;
    auto toVertex = temp.to_vertex_handle(temp.next_halfedge_in_halfface(heh, hfh));
    auto newVertex = temp.split_edge(temp.edge_handle(heh));
    auto heToCollapse = temp.halfedge(newVertex, toVertex);
    auto oppToCollapse_nc = temp.halfedge(newVertex, oppVh_nc);

    if(link_condition(temp, heToCollapse) && heh.is_valid()){
        mesh_.split_edge(mesh_.edge_handle(heh));
        for(auto ch: mesh_.halfedge_cells(oppToCollapse_nc)){
            _cellsAdded->push_back(ch);
        }
        mesh_.collapse_edge(heToCollapse);
//        std::cout << "flip 2-2" << std::endl;
    }
}

// ----------- ** Contraction pass ** ---------------------

void TetLoop::edge_contraction_pass(PriorityQueue* _A){
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

// ----------- ** Insertion pass ** -----------------------

void TetLoop::insertion_pass(PriorityQueue* _A){
//    for each inverted tetrahedron 𝑐 ∈ 𝐶
//    if 𝑐 ∉ 𝑆 for all 𝑆 ∈ G                     𝑐 not yet covered
//    𝑆 ∗ = ⟨{𝑐}⟩                                spawn star at seed 𝑐
//    do                                          growing
//    Choose next tetrahedron 𝑐 ∗ by heuristic
//    𝑆 ∗ = ⟨𝑆 ∗ ∪ {𝑐 ∗ } ⟩
//    if 𝑐 ∗ ∈ 𝑆 for any 𝑆 ∈ G                   ⊲ star collision
//    𝑆 ∗ = ⟨𝑆 ∗ ∪ 𝑆 ⟩                           ⊲ absorb other star
//    G = G \ {𝑆 }
//    while 𝑆 ∗ violates the star conditions
//    G = G ∪ {𝑆 ∗ }
//    return G
    bool saveMesh(true);

    for(auto hfh: mesh_.halffaces()){
        cavityEdge_[hfh] = false;
    }

    std::vector<Star> galaxy;
    while(!_A->empty()){
        auto top = _A->top();
        _A->pop();
        auto ch = top.cell_handle_;
        if(mesh_.is_deleted(ch))
            continue;
        int maxIter = 100;
        bool covered = false;

        for(auto star: galaxy){
            if(std::find(star.tets_.begin(), star.tets_.end(), ch) != star.tets_.end()){
                covered = true;
                break;
            }
        }
        if(covered){
            continue;
        }
        Star newStar(std::set<CellHandle> {ch});
        for(auto hfh: mesh_.cell_halffaces(ch)){
            cavityEdge_[hfh] = true;
        }
        do {
            auto next = findNextCell(newStar);
            if(next.is_valid()){
                newStar.tets_.insert(next);
                for(auto hfh: mesh_.cell_halffaces(next)){
                    cavityEdge_[hfh] = true;
                }
                // check for star collisions
                std::set<int> toRemove;
                for(int i = 0; i < galaxy.size(); ++i){
                    auto star = galaxy[i];
                    // check if the star has already been absorbed
                    if(toRemove.find(i) == toRemove.end() && std::find(star.tets_.begin(), star.tets_.end(), next) != star.tets_.end()){
                        // absorb star
                        newStar.tets_.insert(star.tets_.begin(), star.tets_.end());
                        toRemove.insert(i);
                    }
                }
                for(int id: toRemove){
                    galaxy.erase(galaxy.begin() + id);
                }
//                std::cout << "Boundary: "<< newStar.bounds_.size() << std::endl;
            }
        // star conditions -> centre de chebyshev + injectivity check (* shaped dans 2 domaines)
        } while (checkStarConditions(newStar) && maxIter-- > 0);
//        std::cout << "Boundary: "<< newStar.bounds_.size() << std::endl;
        galaxy.push_back(newStar);

    }

    for(auto star: galaxy){
//        std::cout << "Center: "<< star.center_ << std::endl;
        if(!checkCenter(star)){
            std::cout << "Center of star check failed" << std::endl;
            if(saveMesh){
                createBoundaryMesh(star.bounds_);
            }
            continue;
        }
        for(auto ch: star.tets_){
            if(mesh_.is_deleted(ch)){
                std::cout << "Already deleted cell: " << ch << std::endl;
                continue;
            }
            mesh_.delete_cell(ch);
        }
    }

    // Remplir la galaxy
    // for each star
    // 1. compute center of cheby and add point
    // 2. fill with vertices from boundary (halfface.get_vertices + newPoint)
    for(auto star: galaxy){
        auto newVertex = mesh_.add_vertex(star.center_);
        for(auto bound: star.bounds_){
            // TODO: rajouter la face supprimée
            if(mesh_.is_deleted(bound))
                continue;
            auto vertices = mesh_.get_halfface_vertices(bound);
            auto added = mesh_.add_cell(vertices[0], vertices[2], vertices[1], newVertex);
//            std::cout << "Added cell "<< added << " Valid: "<< added.is_valid() << std::endl;
        }

    }
}

void TetLoop::createBoundaryMesh(std::set<HalfFaceHandle>& _cavityBoundary){
    TriMesh mesh;
    static int id = 0;
    for(auto hfh: _cavityBoundary){
        auto ovm_vertexHandles = mesh_.get_halfface_vertices(hfh);
        std::vector<OpenMesh::SmartVertexHandle> om_vertexHandles;
        for(auto vh: ovm_vertexHandles){
            Point p(mesh_.vertex(vh));
            auto v = mesh.add_vertex(p);
            om_vertexHandles.push_back(v);
        }
        mesh.add_face(om_vertexHandles);
    }
    OpenMesh::IO::write_mesh(mesh, LOGS_MESH + "mesh_dump" + std::to_string(id++) + ".obj");

}

bool TetLoop::checkCenter(Star _star){
    bool pass(false);
    double treshold = 1;
    double dist = 0;
    double averageDist = 10;
    for(auto hfh: _star.bounds_){
        auto face = mesh_.face_handle(hfh);
        auto barycenter = mesh_.barycenter(face);
        dist += (barycenter - _star.center_).norm();
    }
    averageDist = dist / _star.bounds_.size();
//    std::cout << "Average dist: " << averageDist << std::endl;
    pass = averageDist < treshold;

    return pass;

}

bool TetLoop::checkStarConditions(Star _star){
    bool isValid(false);
    _star.bounds_.clear();
    findCavityBoundary(_star.bounds_, _star);
    std::vector<HalfFaceHandle> bounds(_star.bounds_.begin(), _star.bounds_.end());
    isValid = find_chebyshev_center(bounds, CHEBY_THRESHOLD, _star.center_);
//    std::cout << "Center: "<< _star.center_ <<"\nBounds: "<< _star.bounds_.size() << std::endl;

    return isValid;

}

void TetLoop::findCavityBoundary(std::set<HalfFaceHandle>& _cavityBoundary){
    for(auto hfh: mesh_.halffaces()){
        auto opp = mesh_.opposite_halfface_handle(hfh);
        if(cavityEdge_[hfh] && !cavityEdge_[opp] && opp.is_valid()){
            _cavityBoundary.insert(opp);
        }

    }

}

void TetLoop::findCavityBoundary(std::set<HalfFaceHandle>& _cavityBoundary, Star _constraint){
    for(auto tet: _constraint.tets_){
        for(auto ch: mesh_.cell_cells(tet)){
            for(auto hfh: mesh_.cell_halffaces(ch)){
                auto opp = mesh_.opposite_halfface_handle(hfh);
                if(cavityEdge_[hfh] && !cavityEdge_[opp] && opp.is_valid()){
                    _cavityBoundary.insert(opp);
                }
            }
        }
    }
    // This means that only the cavity is composed of one cell only
    if(_cavityBoundary.size() == 0 && _constraint.tets_.size() == 1){
        for(auto hfh: mesh_.cell_halffaces(*_constraint.tets_.begin())){
            auto opp = mesh_.opposite_halfface_handle(hfh);
            // Additional security check
            if(cavityEdge_[hfh] && opp.is_valid()){
                _cavityBoundary.insert(opp);
            }
        }

    }

}

CellHandle TetLoop::findNextCell(Star& _startStar){
    bool printDebug(false);
    CellHandle nextCell(-1);
    _startStar.bounds_.clear();
    findCavityBoundary(_startStar.bounds_, _startStar);
    std::vector<HalfFaceHandle> bounds(_startStar.bounds_.begin(), _startStar.bounds_.end());
    find_chebyshev_center(bounds, CHEBY_THRESHOLD, _startStar.center_);

    int max = -1;
    HalfFaceHandle bestCandidate(-1);
    for(auto hfh: _startStar.bounds_){
        double f_star = -1;
        auto normal = mesh_.normal(hfh);
        auto a = mesh_.barycenter(mesh_.face_handle(hfh));
        f_star = normal.dot(_startStar.center_) - normal.dot(a);
        if(f_star > max){
            bestCandidate = hfh;
            max = f_star;
        }
    }
    if(!bestCandidate.is_valid()){
        if(printDebug)
            std::cout << "Best candidate is not valid" << std::endl;
        return nextCell;
    }

    nextCell = mesh_.incident_cell(bestCandidate);

    if(printDebug){
        if(nextCell.is_valid()){
            std::cout << "Size of boundary: "<< _startStar.bounds_.size() <<"\nNext best cell is "<< nextCell << std::endl;
        }else{
            std::cout << "No valid cell found" << std::endl;
        }

    }

    return nextCell;
}

void TetLoop::triangle_normal_and_centroid(const HalfFaceHandle& hf,
                                 ACG::Vec3d& normal,
                                 ACG::Vec3d& centroid){

    auto face_vertices = mesh_.get_halfface_vertices(hf);
    auto c_evi  = mesh_.vertex(face_vertices[1]) - mesh_.vertex(face_vertices[0]);
    auto c_evi1 = mesh_.vertex(face_vertices[2]) - mesh_.vertex(face_vertices[0]);


    normal = (c_evi.cross(c_evi1));

    auto normal_norm = normal.norm();
    if(normal_norm){
        normal /= normal_norm;
    }

    centroid = (mesh_.vertex(face_vertices[0]) +
                mesh_.vertex(face_vertices[1]) +
                mesh_.vertex(face_vertices[2]))/3.f;

}


bool TetLoop::find_chebyshev_center(const std::vector<HalfFaceHandle>& constraint_hfs,
                          const double radius_lower_bound,
                          ACG::Vec3d& new_position){

    bool print_debug(false);

    new_position = {0,0,0};

    using Program  = CGAL::Quadratic_program<double>;
    using Solution = CGAL::Quadratic_program_solution<double>;

    int max_precision(0);
    int n_constraints(0);

    CGAL::Quadratic_program_options options;
    options.set_pricing_strategy(CGAL::QP_BLAND);     // Bland's rule

    //create solver
    // by default, we have an LP with Ax <= b
    Program lp(CGAL::SMALLER, false, 0, false, 0);

    // now set the non-default entries
    const int X(0), Y(1), Z(2), R(3);

    for(auto hf: constraint_hfs){

            ACG::Vec3d normal, centroid;

            triangle_normal_and_centroid(hf,
                                         normal, centroid);

            auto b = normal.dot(centroid);

            //Ai^T . x + |Ai| . r <= Bi
            lp.set_a(X, n_constraints, normal[0]);
            lp.set_a(Y, n_constraints, normal[1]);
            lp.set_a(Z, n_constraints, normal[2]);

            lp.set_a(R, n_constraints, norm(normal));

            lp.set_b(n_constraints, b);

            n_constraints++;
    }


    if(print_debug){
        std::cout<<" - added "<<n_constraints<<" constraints "<<std::endl;
    }


    //r >= 0
    lp.set_l(R, true, radius_lower_bound);

    // Objective function: -r = (0x + 0y + 0z - 1r)
    lp.set_c(X,  0);
    lp.set_c(Y,  0);
    lp.set_c(Z,  0);
    lp.set_c(R, -1);


    //solve LP

    //std::cout<<" - about to solve system..."<<std::endl;
    Solution s = CGAL::solve_linear_program(lp, double(), options);

    double radius;
    int i(0);
    for(auto val_it = s.variable_values_begin(); val_it != s.variable_values_end(); val_it++){

        double sol_value = (val_it->numerator()/ val_it->denominator());

        if(i<3){
            new_position[i] = sol_value;
        }else{
            radius = sol_value;
        }
        i++;
    }


    if(s.is_infeasible() || s.is_unbounded() || !s.is_valid()){
        if(print_debug){
            std::cout<<" --> infeasible: "<<s.is_infeasible()<<std::endl;
            std::cout<<" -->  unbounded: "<<s.is_unbounded()<<std::endl;
            std::cout<<" -->      valid: "<<s.is_valid()<<std::endl;
        }
        return false;
    }

    return true;

}

// ----------- ** Smoothing pass ** -----------------------

void TetLoop::smoothing_pass(PriorityQueue* _A, int _iterations){
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
