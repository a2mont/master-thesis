#include "TetLoop.hh"

void TetLoop::loop(int _max_iter){
    auto begin = std::chrono::high_resolution_clock::now();
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
            improve_mesh(worstTets);
        }
        reset_queue(quality_queue_);
    }
    if(includeLogs_){
        // register the worst element to stats
        computeQuality();
        addToStats(Stats::TIMESTEP, quality_queue_.top().quality_);
    }
    displayIterationTime(begin, "Complete remeshing");
    if(includeLogs_){
        logStats(stats_, *logger_);
        stats_.flush();
    }

}

void TetLoop::improve_mesh(PriorityQueue& _badTets){
    // B <- set of triangles in M with q < q_min
    // foreach t in B
    // if t still exists and q < q_min improve_tet(M,t,q_min)
    int count(0);
    int total = _badTets.size();
    while (!_badTets.empty()) {
        if(++count % 10 == 0){
            std::cout << "Tet " <<count << "/" << total<< std::endl;
        }
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
    bool topo_active(false), contraction_active(false), insertion_active(false), smoothing_active(false);
    bool changed;
    bool printDebug(false);
    if(printDebug){
        std::cout << "Improve tet: "<< _t.cell_handle_ << std::endl;
    }
    auto begin = std::chrono::high_resolution_clock::now();
    for (int i = 0; i < 1; ++i) {
        int maxIter = 10;
        changed = false;
        // A <- topological_pass(A,M)
        if(topo_active){
            if(printDebug){
                std::cout << "\033[1;36mTopological pass with A of size: \033[0m" << A.size()<< std::endl;
            }
            auto topo_begin = std::chrono::high_resolution_clock::now();
            double qualityBefore(0.), qualityDelta(0.);
            if(includeLogs_){
                computeQuality();
                qualityBefore = quality_queue_.top().quality_;
            }
            do {
                changed = topologial_pass(A);
                if(A.top().quality_ >= q_min_ || A.empty()){
                    if(includeLogs_){
                        qualityDelta = computeQualityDelta(qualityBefore);
                        addToStats(Stats::TOPOLOGY, qualityDelta);
                    }
                    if(printDebug){
                        std::cout << "Exited after topological pass with\n"
                                     "\t-Size of A: "<< A.size() << " -Quality: " <<A.top().quality_<< std::endl;
                        displayIterationTime(topo_begin, "Topological pass");
                    }
                    return;
                }
            } while (changed && maxIter-- > 0);
            if(includeLogs_){
                qualityDelta = computeQualityDelta(qualityBefore);
                addToStats(Stats::TOPOLOGY, qualityDelta);
            }
            if(printDebug){
                displayIterationTime(topo_begin, "Topological pass");
            }
        }
        // A <- edge_contraction_pass(A,M)
        if(contraction_active){
            if(printDebug){
                std::cout << "\033[1;36mEdge contraction pass with A of size: \033[0m" << A.size() << std::endl;
            }
            auto contraction_begin = std::chrono::high_resolution_clock::now();
            double qualityBefore(0.), qualityDelta(0.);
            if(includeLogs_){
                computeQuality();
                qualityBefore = quality_queue_.top().quality_;
            }
            edge_contraction_pass(A);
            if(A.top().quality_ >= q_min_ || A.empty()){
                if(includeLogs_){
                    qualityDelta = computeQualityDelta(qualityBefore);
                    addToStats(Stats::CONTRACTION, qualityDelta);
                }
                if(printDebug){
                    std::cout << "Exited after contraction pass with\n"
                                 "\t-Size of A: "<< A.size() << " -Quality: " <<A.top().quality_<< std::endl;
                    displayIterationTime(contraction_begin, "Contraction pass");
                }
                return;
            }
            if(includeLogs_){
                qualityDelta = computeQualityDelta(qualityBefore);
                addToStats(Stats::CONTRACTION, qualityDelta);
            }
            if(printDebug){
                displayIterationTime(contraction_begin, "Contraction pass");
            }

        }
        // A <- insertion_pass(A,M)
        if(insertion_active){
            if(printDebug){
                std::cout << "\033[1;36mInsertion pass with A of size: \033[0m" << A.size()<< std::endl;
            }
            auto insertion_begin = std::chrono::high_resolution_clock::now();
            double qualityBefore(0.), qualityDelta(0.);
            if(includeLogs_){
                computeQuality();
                qualityBefore = quality_queue_.top().quality_;
            }
            insertion_pass(A);
            if(A.top().quality_ >= q_min_ || A.empty()){
                if(includeLogs_){
                    qualityDelta = computeQualityDelta(qualityBefore);
                    addToStats(Stats::INSERTION, qualityDelta);
                }
                if(printDebug){
                    std::cout << "Exited after insertion pass with\n"
                                 "\t-Size of A: "<< A.size() << " -Quality: " <<A.top().quality_<< std::endl;
                    displayIterationTime(insertion_begin, "Insertion pass");
                }
                return;
            }
            if(includeLogs_){
                qualityDelta = computeQualityDelta(qualityBefore);
                addToStats(Stats::INSERTION, qualityDelta);
            }
            if(printDebug){
                displayIterationTime(insertion_begin, "Insertion pass");
            }
        }
        // smoothing_pass
        if(smoothing_active){
            if(printDebug){
                std::cout << "\033[1;36mSmoothing pass with A of size: \033[0m" << A.size()<< std::endl;
            }
            auto smoothing_begin = std::chrono::high_resolution_clock::now();
            double qualityBefore(0.), qualityDelta(0.);
            if(includeLogs_){
                computeQuality();
                qualityBefore = quality_queue_.top().quality_;
            }
            smoothing_pass(A);
            if(A.top().quality_ >= q_min_ || A.empty()){
                if(includeLogs_){
                    qualityDelta = computeQualityDelta(qualityBefore);
                    addToStats(Stats::SMOOTHING, qualityDelta);
                }
                if(printDebug){
                    std::cout << "Exited after smoothing pass with\n"
                                 "-A: "<< A.size() << " , Quality: " <<A.top().quality_<< std::endl;
                    displayIterationTime(smoothing_begin, "Smoothing pass");
                }
                return;
            }
            if(includeLogs_){
                qualityDelta = computeQualityDelta(qualityBefore);
                addToStats(Stats::SMOOTHING, qualityDelta);
            }
            if(printDebug){
                displayIterationTime(smoothing_begin, "Smoothing pass");
            }
        }
    }
    if(printDebug){
        displayIterationTime(begin, "Tet " + std::to_string(_t.cell_handle_.idx()));
    }

}

// ----------- ** Topological pass ** ---------------------

bool TetLoop::topologial_pass(PriorityQueue& _A){
    //flips,edge removal, multi-face removal operations
    bool changed, edge, face;
    changed = edge = face = false;
    std::vector<CellHandle> cellsAdded;
    std::vector<int> counter({0,0,0});
    int total_face = 0;
    int total_edge = 0;
    // for t in A that still exists
    //      for each edge e of t (if t still exists)
    //          attempt to remove e
    //      for each f of t (if t still exists)
    //          attempt to remove f (multi-face, 2-3/2-2 flip)
    // A <- A U new tets created by the operation
    while(!_A.empty()){
        auto top = _A.top();
        _A.pop();
        auto ch = top.cell_handle_;
        //        std::cout << top.quality_<< " vs " <<QualityEvaluation::evaluate(ch,mesh_) << std::endl;
        if(mesh_.is_deleted(ch)) continue;

        std::vector<EdgeHandle> edges;
        std::vector<FaceHandle> faces;
        for(auto c_eh: mesh_.cell_edges(ch)){
            edges.push_back(c_eh);
        }
        auto edge_begin = std::chrono::high_resolution_clock::now();
        for(auto e: edges){
            if(!mesh_.is_deleted(e)){
                edge = edgeRemoval(e, cellsAdded,true);
                changed = changed || edge;
                total_edge = edge ? ++total_edge : total_edge;
            }
        }
        //        displayIterationTime(edge_begin, "Edge part");
        if(changed || mesh_.is_deleted(ch)) continue;

        for(auto c_fh: mesh_.cell_faces(ch)){
            faces.push_back(c_fh);
        }
        auto face_begin = std::chrono::high_resolution_clock::now();
        for(auto f: faces){
            if(!mesh_.is_deleted(f)){
                face = faceRemoval(f, cellsAdded,counter, false);
                changed = changed || face;
                total_face = face ? ++total_face : total_face;
            }
        }
        //        displayIterationTime(face_begin, "Face part");
        // if nothing happened, we push the cell to the queue again
        if(!changed && !mesh_.is_deleted(ch)){
            cellsAdded.push_back(ch);
        }
    }
    std::cout << "Total changes edges: " << total_edge << std::endl;

    std::cout << "Method used: \n\t2-3: "<< counter[0] <<
                 "\t2-2: "<< counter[1] << "\tMulti: "<< counter[2]<< "\nout of "<< total_face << " changes "
              << std::endl;

    for(auto ch: cellsAdded){
        if(mesh_.is_deleted(ch)){
            continue;
        }
        double quality = QualityEvaluation::evaluate(ch, mesh_);
        Tet toAdd(ch, quality);
        _A.push(toAdd);
    }
    return changed;
}

// Call when we don't track the cells added by the process (e.g. Tests.cc)
bool TetLoop::edgeRemoval(EdgeHandle _eh, bool _verbose){
    std::vector<CellHandle> ignore;
    return edgeRemoval(_eh, ignore, _verbose);
}

bool TetLoop::edgeRemoval(EdgeHandle _eh, std::vector<CellHandle>& _cellsAdded, bool _verbose){
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

    // find the new hafledges between new vertex and one ring
    for(auto vheh: tempBaseMesh.outgoing_halfedges(newVertex)){
        if(tempBaseMesh.to_vertex_handle(vheh) == to || tempBaseMesh.to_vertex_handle(vheh) == from)
            continue;
        newHalfEdges.push_back(vheh);
    }
    // find the best halfedge to collapse
    for(auto heh: newHalfEdges){
        auto tempMesh = tempBaseMesh;
        if(!link_condition(tempMesh, heh)) continue;
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

    // collapse the best halfedge on the real mesh
    if(maxQuality > q_old){
        if(link_condition(mesh_, toCollapse) && toCollapse.is_valid()){
            if(useWorldSpace_){
                auto wm_copy = world_mesh_;
                wm_copy.split_edge(_eh);
                wm_copy.collapse_edge(toCollapse);
                bool valid = validateWorldMesh(wm_copy);
                // operation is not valid in world space -> abort
                if(!valid){
                    if(_verbose){
                        std::cout << "Edge removal aborted in world space!!!" << std::endl;
                    }
                    return changed;
                }
                // update world space
                world_mesh_.split_edge(_eh);
                auto newFrom_wm = world_mesh_.collapse_edge(toCollapse);
                Smoothing::smooth(world_mesh_, newFrom_wm);
            }
            std::vector<CellHandle> adjacentCells;
            auto addedVertex = mesh_.split_edge(_eh);
            for(auto vch: mesh_.vertex_cells(addedVertex)){
                adjacentCells.push_back(vch);
            }
            auto newFrom = mesh_.collapse_edge(toCollapse);
            Smoothing::smooth(mesh_, newFrom);
            // find the altered tets
            for(auto ch: adjacentCells){
                if(mesh_.is_deleted(ch)) continue;
                _cellsAdded.push_back(ch);
            }
            changed = true;
            //            std::cout << "Tets added: "<< _cellsAdded.size() << " , changed:"<< changed << std::endl;
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

// Call when we don't track the cells added by the process (e.g. Tests.cc)
bool TetLoop::faceRemoval(FaceHandle _fh, bool _verbose){
    std::vector<CellHandle> ignore;
    std::vector<int> ig;
    return faceRemoval(_fh, ignore, ig, _verbose);
}

bool TetLoop::faceRemoval(FaceHandle _fh, std::vector<CellHandle>& _cellsAdded, std::vector<int>& _counter,
                          bool _verbose){
    bool changed(false);
    auto tempMesh = mesh_;
    PriorityQueue tempQueue;

    /* track the changes of the different remeshings
     * ids:
     *  - 0: 2-3 split
     *  - 1: multiface
     *  - 2: Actual result
    */
    std::vector<std::vector<CellHandle>> addedTets(3);
    int bestId = -1;

    // Quality before any operation
    computeQuality(tempQueue, mesh_);
    double q_old = tempQueue.top().quality_;
    double q_max = q_old;
    double q_new = 0;

    flip23(tempMesh,_fh, addedTets[0]);
    // Quality after 2-3 flip
    computeQuality(tempQueue, tempMesh);
    q_new = tempQueue.top().quality_;
    std::cout << "2-3\nNew: " << q_new << " Max: "<< q_max << std::endl;
    if(q_new > q_max){
        q_max = q_new;
        changed = true;
        bestId = 0;
        if(_verbose){
            std::cout << "2-3 flip has better quality\n" <<
                         "added tets: "<< addedTets[0].size() << std::endl;
        }
    }
    tempMesh = mesh_;

    multiFace(tempMesh,_fh, addedTets[1]);
    computeQuality(tempQueue, tempMesh);
    q_new = tempQueue.top().quality_;
    std::cout << "Multi\nNew: " << q_new << " Max: "<< q_max << std::endl;
    if(q_new > q_max){
        q_max = q_new;
        changed = true;
        bestId = 1;
        if(_verbose){
            std::cout << "Multi face has better quality\n" <<
                         "Added tets: "<< addedTets[1].size() << std::endl;
        }
    }
    if(bestId == -1){
        return changed = false;
    }
    ++_counter[bestId];

    switch(bestId){
    case 0:
        flip23(mesh_, _fh, addedTets[2]);
        break;
    case 1:
        multiFace(mesh_, _fh, addedTets[2]);
        break;

    }

    for(auto ch: addedTets[2]){
        if(!mesh_.is_deleted(ch)){
            _cellsAdded.push_back(ch);
        }
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
void TetLoop::multiFace(TetrahedralMesh& _mesh, FaceHandle _fh, std::vector<CellHandle>& _cellsAdded){
    if(_mesh.is_boundary(_fh))
        return;
    auto hfhs = _mesh.face_halffaces(_fh);
    auto a = _mesh.halfface_opposite_vertex(hfhs[0]);
    auto b = _mesh.halfface_opposite_vertex(hfhs[1]);
    std::vector<VertexHandle> faceVertices; // u,v,w

    for(auto vh: _mesh.get_halfface_vertices(hfhs[0])){
        faceVertices.push_back(vh);
    }
    auto results_uv = testNeighbor(_mesh,a,b,faceVertices[0],faceVertices[1]);
    auto results_vw = testNeighbor(_mesh,a,b,faceVertices[1],faceVertices[2]);
    auto results_wu = testNeighbor(_mesh,a,b,faceVertices[2],faceVertices[0]);
    double q_old = std::min({QualityEvaluation::evaluate(a,faceVertices[0],faceVertices[1],faceVertices[2], _mesh),
                             QualityEvaluation::evaluate(faceVertices[0],faceVertices[1],faceVertices[2],b, _mesh),
                             results_uv.o_,
                             results_vw.o_,
                             results_wu.o_
                            });
    double q_new = std::min({results_uv.n_, results_vw.n_, results_wu.n_});
    if(q_new > q_old){
        flip23(_mesh, _fh, _cellsAdded);
        std::vector<FaceWithChildren> total;
        total.insert(total.end(), results_uv.h_.begin(), results_uv.h_.end());
        total.insert(total.end(), results_vw.h_.begin(), results_vw.h_.end());
        total.insert(total.end(), results_wu.h_.begin(), results_wu.h_.end());
        for(auto& g: total){
            flip32Recurse(_mesh, g, _fh, _cellsAdded);
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
TetLoop::TestNeighborResult TetLoop::testNeighbor(TetrahedralMesh& _mesh,
                                                  VertexHandle a,
                                                  VertexHandle b,
                                                  VertexHandle u,
                                                  VertexHandle w){
    double q_uw = QualityEvaluation::evaluate(a,b,u,w,_mesh);
    auto he_uw = _mesh.find_halfedge(u,w);
    auto uw = _mesh.edge_handle(he_uw);
    FaceWithChildren g;
    VertexHandle v;
    bool correctBoundary = !_mesh.is_boundary(uw) ||
            (_mesh.is_boundary(a) && _mesh.is_boundary(b));

    if(_mesh.valence(uw) == 4 && correctBoundary){
        for(auto hehfh: _mesh.halfedge_halffaces(he_uw)){
            auto next_he = _mesh.next_halfedge_in_halfface(he_uw, hehfh);
            if(_mesh.to_vertex_handle(next_he) == a || _mesh.to_vertex_handle(next_he) == b){
                continue;
            }
            g = _mesh.face_handle(hehfh);
            v = _mesh.to_vertex_handle(next_he);
            break;
        }
        double j_uv = orient3D(_mesh,a,b,u,v);
        double j_vw = orient3D(_mesh,a,b,v,w);
        double j_wu = orient3D(_mesh,a,b,w,u);
        if((j_uv > 0 && j_vw > 0) ||
                (j_vw > 0 && j_wu > 0) ||
                (j_wu > 0 && j_uv > 0)){
            auto result_uv = testNeighbor(_mesh,a,b,u,v);
            auto result_vw = testNeighbor(_mesh,a,b,v,w);
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
            double q_old = std::min({QualityEvaluation::evaluate(a,u,v,w,_mesh),
                                     QualityEvaluation::evaluate(u,v,w,b,_mesh),
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

double TetLoop::orient3D(TetrahedralMesh& _mesh,
                         VertexHandle _a,
                         VertexHandle _b,
                         VertexHandle _c,
                         VertexHandle _d){
    std::vector<Point> points = {_mesh.vertex(_a), _mesh.vertex(_b), _mesh.vertex(_c), _mesh.vertex(_d)};
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

void TetLoop::flip32Recurse(TetrahedralMesh& _mesh,
                            TetLoop::FaceWithChildren _g,
                            FaceHandle _parent,
                            std::vector<CellHandle>& _cellsAdded){
    EdgeHandle toRemove;
    // find common edge with parent
    for(auto parentEdge: _mesh.face_edges(_parent)){
        for(auto edge: _mesh.face_edges(_g.fh_)){
            if(edge == parentEdge){
                toRemove = edge;
            }
        }
    }
    flip32(_mesh, toRemove, _cellsAdded);
    for(auto h: _g.children_){
        flip32Recurse(_mesh, h, _g.fh_, _cellsAdded);
    }
}

void TetLoop::flip32(TetrahedralMesh& _mesh,
                     EdgeHandle _eh,
                     std::vector<CellHandle>& _cellsAdded){
    if(_mesh.is_deleted(_eh)) return;
    bool printDebug(false);
    auto ends = _mesh.edge_vertices(_eh);
    HalfEdgeHandle toCollapse(-1);
    auto temp = _mesh;
    VertexHandle newVertex(-1);
    try {
        newVertex = temp.split_edge(_eh);
    } catch (...) {
        std::cout << "------------ Error splitting edge: "<< _eh << "--------------" << std::endl;
    }
    if(!newVertex.is_valid() || _mesh.is_boundary(ends[0]) || _mesh.is_boundary(ends[1])){
        return;
    }
    for(auto heh: _mesh.outgoing_halfedges(newVertex)){
        if(_mesh.to_vertex_handle(heh) != ends[0] && _mesh.to_vertex_handle(heh) != ends[1]){
            toCollapse = heh;
        }
    }
    if(link_condition(temp,toCollapse)){
        if(useWorldSpace_){
            auto wm_copy = world_mesh_;
            wm_copy.split_edge(_eh);
            wm_copy.collapse_edge(toCollapse);
            bool valid = validateWorldMesh(wm_copy);
            // operation is not valid in world space -> abort
            if(!valid){
                if(printDebug){
                    std::cout << "3-2 split aborted in world space!!!" << std::endl;
                }
                return;
            }
            // update world space
            world_mesh_.split_edge(_eh);
            world_mesh_.collapse_edge(toCollapse);
        }
        auto addedVertex = _mesh.split_edge(_eh);
        std::set<CellHandle> adjacentCells;
        for(auto ch: _mesh.halfedge_cells(toCollapse)){
            adjacentCells.insert(ch);
        }
        for(auto ch: _mesh.vertex_cells(addedVertex)){
            // checks that the cells aren't adjacent to the collapsed he
            if(adjacentCells.find(ch) == adjacentCells.end()){
                _cellsAdded.push_back(ch);
            }
        }
        _mesh.collapse_edge(toCollapse);
    }
}

void TetLoop::flip23(TetrahedralMesh& _mesh,
                     FaceHandle _fh,
                     std::vector<CellHandle>& _cellsAdded){
    if(_mesh.is_deleted(_fh)) return;
    bool printDebug(false);
    auto hfh = _mesh.face_halffaces(_fh)[0];
    auto toVertex = _mesh.halfface_opposite_vertex(hfh);
    if(_mesh.is_boundary(toVertex)){
        hfh = _mesh.face_halffaces(_fh)[1];
        toVertex = _mesh.halfface_opposite_vertex(hfh);
    }
    if(_mesh.is_boundary(toVertex)) return;

    auto hfh_opp = _mesh.opposite_halfface_handle(hfh);
    auto vh_opp = _mesh.halfface_opposite_vertex(hfh_opp);

    auto temp = _mesh;
    auto newVertex = temp.split_face(_fh);
    auto heToCollapse = temp.find_halfedge(newVertex, toVertex);
    auto remainingHe = temp.find_halfedge(vh_opp, newVertex);

    if(link_condition(temp, heToCollapse)){
        if(useWorldSpace_){
            auto wm_copy = world_mesh_;
            wm_copy.split_face(_fh);
            wm_copy.collapse_edge(heToCollapse);
            bool valid = validateWorldMesh(wm_copy);
            // operation is not valid in world space -> abort
            if(!valid){
                if(printDebug){
                    std::cout << "2-3 split aborted in world space!!!" << std::endl;
                }
                return;
            }
            // update world space
            world_mesh_.split_face(_fh);
            world_mesh_.collapse_edge(heToCollapse);
        }
        _mesh.split_face(_fh);
        for(auto ch: _mesh.halfedge_cells(remainingHe)){
            _cellsAdded.push_back(ch);
        }
        _mesh.collapse_edge(heToCollapse);

    }
}

void TetLoop::flip22(TetrahedralMesh& _mesh,
                     FaceHandle _fh,
                     std::vector<CellHandle>& _cellsAdded){
    if(_mesh.is_deleted(_fh)) return;
    bool printDebug(false);
    HalfEdgeHandle heh(-1);
    HalfEdgeHandle toCollapse(-1);
    auto hfh = _mesh.face_halffaces(_fh)[0];
    for(auto he: _mesh.face_halfedges(_fh)){
        if(!_mesh.is_boundary(he)){
            heh = he;
            break;
        }
    }
    // Get the heh going in the opposite direction from collapse (_nc = next cell)
    auto adjHf = _mesh.adjacent_halfface_in_cell(hfh, heh);
    auto oppHf_nc = _mesh.opposite_halfface_handle(adjHf);
    auto oppHe_nc = _mesh.opposite_halfedge_handle(heh);
    if(!oppHf_nc.is_valid() || !oppHe_nc.is_valid()){
        return;
    }
    auto adjHf_nc = _mesh.adjacent_halfface_in_cell(oppHf_nc, oppHe_nc);
    auto next_He = _mesh.next_halfedge_in_halfface(oppHe_nc, adjHf_nc);

    auto oppVh_nc = _mesh.to_vertex_handle(next_He);
    auto temp = _mesh;
    auto toVertex = temp.to_vertex_handle(temp.next_halfedge_in_halfface(heh, hfh));
    auto newVertex = temp.split_edge(temp.edge_handle(heh));
    auto heToCollapse = temp.find_halfedge(newVertex, toVertex);
    auto oppToCollapse_nc = temp.find_halfedge(newVertex, oppVh_nc);

    if(link_condition(temp, heToCollapse) && heh.is_valid()){
        auto ehToSplit = _mesh.edge_handle(heh);
        if(useWorldSpace_){
            auto wm_copy = world_mesh_;
            wm_copy.split_edge(ehToSplit);
            wm_copy.collapse_edge(heToCollapse);
            bool valid = validateWorldMesh(wm_copy);
            // operation is not valid in world space -> abort
            if(!valid){
                if(printDebug){
                    std::cout << "2-2 split aborted in world space!!!" << std::endl;
                }
                return;
            }
            // update world space
            world_mesh_.split_edge(ehToSplit);
            world_mesh_.collapse_edge(heToCollapse);
        }

        _mesh.split_edge(ehToSplit);
        for(auto ch: _mesh.halfedge_cells(oppToCollapse_nc)){
            _cellsAdded.push_back(ch);
        }
        _mesh.collapse_edge(heToCollapse);
    }
}

// ----------- ** Contraction pass ** ---------------------

void TetLoop::edge_contraction_pass(PriorityQueue& _A){
    bool printDebug(false);
    std::vector<EdgeHandle> edges;
    std::vector<CellHandle> newTets;
    PriorityQueue tempA = _A;
    TetrahedralMesh mesh_copy = mesh_;
    computeQuality();
    double q_old = quality_queue_.top().quality_;
    while(!_A.empty()){
        bool changed(false);
        edges.clear();
        auto top = _A.top();
        _A.pop();
        auto ch = top.cell_handle_;
        // E <- edges of the tet in A
        for(auto eh: mesh_.cell_edges(ch)){
            edges.push_back(eh);
        }
        for(auto e: edges){
            if(mesh_.is_deleted(e)){
                continue;
            }
            auto remain = contractEdge(e, newTets);
            changed = changed || remain.is_valid();
        }
        // if nothing happenend to the tet, push back onto queue
        if(!changed){
            newTets.push_back(ch);
        }

    }
    if(newTets.empty()){
        cleanQualityQueue(tempA, mesh_);
        _A = tempA;
        if(printDebug){
            std::cout << "new cells empty" << std::endl;
        }
        return;
    }
    computeQuality();
    double q_new = quality_queue_.top().quality_;
//    if(q_new < q_old){
//        std::cout << "Contraction rolled back q_old: "<< q_old<< " vs q_new: "<< q_new << std::endl;
//        cleanQualityQueue(tempA, mesh_);
//        _A = tempA;
//        mesh_ = mesh_copy; // TODO find some way to make it work
//        return;
//    }
    std::set<CellHandle> newTets_no_dupl;
    // return surviving set in A and the tets in M altered by the contractions
    for(auto new_fh: newTets){
        if(mesh_.is_deleted(new_fh)) continue;
        //remove duplicate from newTets
        auto result = newTets_no_dupl.insert(new_fh);
        if(result.second) {
            _A.push(Tet(new_fh, QualityEvaluation::evaluate(new_fh,mesh_)));
        }
    }
    if(printDebug){
        std::cout << "Total: "<< newTets.size() << " vs duplicates/deleted removed: "<< newTets_no_dupl.size() << std::endl;
    }
}

VertexHandle TetLoop::contractEdge(EdgeHandle _eh, std::vector<CellHandle>& _tetsAltered){
    bool printDebug(false);
    VertexHandle remain = mesh_.InvalidVertexHandle;

    if(mesh_.is_deleted(_eh)){
        return remain;
    }
    auto he = mesh_.edge_halfedges(_eh)[0];
    auto from = mesh_.from_vertex_handle(he);
    auto to = mesh_.to_vertex_handle(he);
    if((mesh_.is_boundary(from) && mesh_.is_boundary(to))){
        return remain;
    }
    auto temp = mesh_;
    // if only the from vertex is on the boundary, invert the collapse direction
    if(mesh_.is_boundary(from)){
        he = mesh_.edge_halfedges(_eh)[1];
        from = mesh_.from_vertex_handle(he);
        to = mesh_.to_vertex_handle(he);
    }
    if(printDebug){
        std::cout << "Collapse vertices: " << from.idx() << " "  << to.idx() << std::endl;
    }

    // used to manually move vertices in testing
    if(constraint_vhs_.count(from.idx()) > 0){
        constraint_vhs_.erase(from.idx());
        constraint_vhs_[to.idx()] = to.idx();
    }
    if(link_condition(mesh_,he)){
        auto tempRemain = temp.collapse_edge(he);
        Smoothing::smooth(temp, tempRemain);
        PriorityQueue tempQueue;
        computeQuality(tempQueue, temp);
        if(tempQueue.top().quality_ > 1e-2){
            if(useWorldSpace_){
                auto wm_copy = world_mesh_;
                wm_copy.collapse_edge(he);
                bool valid = validateWorldMesh(wm_copy);
                // operation is not valid in world space -> abort
                if(!valid){
                    if(printDebug){
                        std::cout << "Edge contraction aborted in world space!!!" << std::endl;
                    }
                    return remain;
                }
                // update world space
                world_mesh_.collapse_edge(he);
            }
            remain = mesh_.collapse_edge(he);
            if(!remain.is_valid()){
                std::cout << "Remain non valid" << std::endl;
                exit(EXIT_FAILURE);
                return remain;
            }
            Smoothing::smooth(mesh_, remain);
            // get the tets around new point
            for(auto v_ch: mesh_.vertex_cells(remain)){
                _tetsAltered.push_back(v_ch);
            }
        }
    }

    return remain;
}

// ----------- ** Insertion pass ** -----------------------


void TetLoop::insertion_pass(PriorityQueue& _A){
    bool saveMesh(false);
    bool printDebug(false);

    for(auto hfh: mesh_.halffaces()){
        cavityEdge_[hfh] = false;
    }
    std::vector<CellHandle> cellsToAdd;
    std::vector<Star> galaxy;
    std::vector<Star> galaxy_wm;
    CellHandle lastAdded(-1);

    auto copy = mesh_;
    auto wm_copy = world_mesh_;
    auto wm_temp = world_mesh_;
    PriorityQueue tempA = _A;
    PriorityQueue queueCopy;
    computeQuality(queueCopy, copy);
    //    for each inverted tetrahedron 𝑐 ∈ 𝐶
    while(!_A.empty()){
        auto top = _A.top();
        _A.pop();
        auto ch = top.cell_handle_;
        if(mesh_.is_deleted(ch)) continue;
        cellsToAdd.push_back(ch);

        int maxStarSize = 5;
        bool covered = false;
        //    if 𝑐 ∉ 𝑆 for all 𝑆 ∈ G                     𝑐 not yet covered
        for(auto& star: galaxy){
            if(std::find(star.tets_.begin(), star.tets_.end(), ch) != star.tets_.end()){
                covered = true;
                break;
            }
        }
        if(covered){
            //            std::cout << "Already covered" << std::endl;
            continue;
        }
        //    𝑆 ∗ = ⟨{𝑐}⟩                                spawn star at seed 𝑐
        Star newStar(std::set<CellHandle> {ch});
        lastAdded = ch;
        for(auto hfh: mesh_.cell_halffaces(ch)){
            auto opp = mesh_.opposite_halfface_handle(hfh);
            cavityEdge_[hfh] = true;
            if(useWorldSpace_){
                cavityEdge_wm_[hfh] = true;
            }
            if(opp.is_valid()){
                newStar.bounds_.insert(opp);
            }
        }
        if(newStar.bounds_.size() == 4){
            find_chebyshev_center(mesh_, newStar.bounds_, CHEBY_THRESHOLD, newStar.center_);
        }else{
            std::cout <<
                         "Error on initial cavity boundary size !!\n\tExpected: 4  Value: "
                      << newStar.bounds_.size()<< std::endl;
        }
        do {
            //    Choose next tetrahedron 𝑐∗ by heuristic
            auto next = findNextCell(newStar, galaxy);
            if(next.is_valid()){
                //    𝑆 ∗ = ⟨𝑆 ∗ ∪ {𝑐 ∗ } ⟩
                newStar.tets_.insert(next);
                lastAdded = next;
                for(auto hfh: mesh_.cell_halffaces(next)){
                    cavityEdge_[hfh] = true;
                    if(useWorldSpace_){
                        cavityEdge_wm_[hfh] = true;
                    }
                }
            }
            // star conditions -> centre de chebyshev + injectivity check (* shaped dans 2 domaines)
        } while (checkStarConditions(newStar, lastAdded) && maxStarSize-- > 0);
        //    G = G ∪ {𝑆 ∗ }
        galaxy.push_back(newStar);
        if(useWorldSpace_){
            galaxy_wm.push_back(newStar);
        }
    }

    for(int i = 0; i < galaxy.size();++i){
        std::cout<<" -------v----------------------------"<<std::endl;
        auto& star = galaxy[i];
        if(!checkCenter(star)){
            if(printDebug){
                std::cout << "Center of star check failed" << std::endl;
            }
            if(saveMesh){
                cavityMesh3D(star);
            }
            // If center fails, cancel
            _A = tempA;
            return;
        }
        for(auto hfh: star.bounds_){
            if(mesh_.is_deleted(hfh)){
                std::cout<<" star bounds hf is deleted"<<std::endl;
            }
            ACG::Vec3d normal, centroid;
            triangle_normal_and_centroid(mesh_, hfh,
                                         normal, centroid);

            if(mesh_.is_boundary(hfh)){
                std::cout<<" added boundary hf "<<vectorToString(mesh_.get_halfface_vertices(hfh))
                        <<" with centroid "<<centroid<<" and normal "<<normal<<std::endl;
            }
            // add vertices for cell reconstruction
            auto vertices = mesh_.get_halfface_vertices(hfh);
            star.reconstructionVectors_.push_back(vertices);

        }
        for(auto ch: star.tets_){
            if(mesh_.is_deleted(ch)){
                std::cout << "Already deleted cell: " << ch << std::endl;
                continue;
            }
            //            if(mesh_.is_boundary(ch)){
            //                for(auto c_hf: mesh_.cell_halffaces(ch)){
            //                    auto opp = mesh_.opposite_halfface_handle(c_hf);
            //                    if(mesh_.is_boundary(opp)){
            //                        auto verts = mesh_.get_halfface_vertices(opp);
            //                        star.surfaceDict_[opp] = std::vector<VertexHandle>(verts);

            //                    }
            //                }
            //            }

            if(useWorldSpace_){
                wm_temp.delete_cell(ch);
                cleanMesh(wm_temp);
                bool valid = validateWorldMesh(wm_temp);
                // operation is not valid in world space -> abort
                if(!valid){
                    if(printDebug){
                        std::cout << "Delete step of insertion aborted in world space!!!" << std::endl;
                    }
                    continue;
                }
            }
            mesh_.delete_cell(ch);
        }
    }

    //    cleanMesh(mesh_, true);

    // Remplir la galaxy
    // for each star
    // 1. compute center of cheby and add point
    // 2. fill with vertices from boundary (halfface.get_vertices + newPoint)
    for(int i = 0; i < galaxy.size(); ++i){
        auto& star = galaxy[i];
        // copy to revert to in case world space aborts operation
        auto tempMesh = mesh_;
        auto newVertex = mesh_.add_vertex(star.center_);
        if(printDebug){
            std::cout << "Added vertex: "
                      << newVertex << " at "<< star.center_
                      <<" (valid = "<< newVertex.is_valid() << ")" << std::endl;
        }
        for(auto vertices: star.reconstructionVectors_){
            if(mesh_.is_deleted(vertices[0])
                    || mesh_.is_deleted(vertices[1])
                    || mesh_.is_deleted(vertices[2])){
                std::cout << "DELETED" << std::endl;
                continue;
            }
            auto added = mesh_.add_cell(vertices[0], vertices[2], vertices[1], newVertex);
            for(auto chf_it = mesh_.chf_iter(added); chf_it.valid(); chf_it++){
                if(mesh_.is_boundary(*chf_it)){
                    ACG::Vec3d normal, centroid;

                    std::cout<<" reconstructed boundary hf "<<vectorToString(mesh_.get_halfface_vertices(*chf_it))
                            <<" with centroid "<<centroid<<" and normal "<<normal<<std::endl;
                }
            }
            Smoothing::smooth(mesh_, newVertex);
            //            std::cout << "Added cell "<< added.idx() << " valid = "<< added.is_valid() << std::endl;
            if(added.is_valid()){
                cellsToAdd.push_back(added);
            }
        }

        //        std::vector<VertexHandle> vertices;
        //        for(auto bound: star.bounds_){
        //            if(mesh_.is_deleted(bound)){
        //                if(star.surfaceDict_.find(bound) == star.surfaceDict_.end()){
        //                    continue;
        //                }
        //                vertices = star.surfaceDict_[bound];

        //            }else{
        //                vertices = mesh_.get_halfface_vertices(bound);
        //            }
        //            // ensure operation is valid in world space
        ////            if(useWorldSpace_){
        ////                Point center_wm = {0,0,0};
        ////                find_chebyshev_center(wm_temp, star.bounds_, CHEBY_THRESHOLD, center_wm);
        ////                auto newVertex_wm = wm_temp.add_vertex(center_wm);
        ////                wm_temp.add_cell(vertices[0], vertices[2], vertices[1], newVertex_wm);
        ////                Smoothing::smooth(wm_temp, newVertex_wm);
        ////                bool valid = validateWorldMesh(wm_temp);
        ////                if(!valid){
        ////                    if(printDebug){
        ////                        std::cout << "Insertion aborted in world space" << std::endl;
        ////                    }
        ////                    // removes the added vertex
        ////                    mesh_ = tempMesh;
        ////                    continue;
        ////                }
        ////                // update world space
        ////                world_mesh_ = wm_temp;
        ////            }
        //            // this will not be executed if world mesh aborts operation
        //            // order = 0,2,1 since boundary is facing opposite to the cavity
        //            auto added = mesh_.add_cell(vertices[0], vertices[2], vertices[1], newVertex);
        //            Smoothing::smooth(mesh_, newVertex);
        //            if(added.is_valid()){
        //                cellsToAdd.push_back(added);
        //            }

        //        }

    }
    // Ensure that no unlinked elements remain
    cleanMesh(mesh_);

    //    computeQuality();
    //    if(quality_queue_.top().quality_ < queueCopy.top().quality_){
    //        mesh_ = copy;
    //        cleanQualityQueue(tempA, mesh_);
    //        _A = tempA;
    //        world_mesh_ = wm_temp;
    //        if(printDebug){
    //            std::cout << "Quality did not improve, revert ("<<
    //                         quality_queue_.top().quality_ << " vs "<<
    //                         queueCopy.top().quality_ << ")"<< std::endl;
    //        }
    //        return;
    //    }

    for(auto ch: cellsToAdd){
        if(mesh_.is_deleted(ch)){
            continue;
        }
        double quality = QualityEvaluation::evaluate(ch, mesh_);
        Tet toAdd(ch, quality);
        _A.push(toAdd);
    }
}


void TetLoop::recoverBadStar(Star& _star, std::string _meshName){
    bool printDebug(false);
    TetrahedralMesh temp;
    OpenVolumeMesh::IO::FileManager fm;
    fm.readFile(LOGS_MESH + _meshName, temp);
    if(temp.n_cells() == 0){
        if(printDebug)
            std::cout << "Error loading mesh !" << std::endl;
        return;
    }
    std::set<HalfFaceHandle> hfs;
    for(auto hf:temp.halffaces()){
        if(temp.is_boundary(hf)){
            hfs.insert(hf);
        }
    }
    ACG::Vec3d new_pos;
    bool success = find_chebyshev_center(temp, hfs, TetLoop::CHEBY_THRESHOLD, new_pos);

    if(success){
        std::cout << "Center moved from "<< _star.center_ << " to " << new_pos << std::endl;
        _star.center_ = new_pos;
        return;
    }
    if(printDebug){
        std::cout << "Could not find cheby center" << std::endl;
    }


}

void TetLoop::compareStarWithMesh(Star& _star, TetrahedralMesh& _mesh){
    findCavityBoundary(_star);
    find_chebyshev_center(mesh_, _star.bounds_, CHEBY_THRESHOLD, _star.center_, true);

    std::set<HalfFaceHandle> hfs;
    for(auto hf:_mesh.halffaces()){
        if(_mesh.is_boundary(hf)){
            hfs.insert(hf);
        }
    }
    ACG::Vec3d new_pos;
    std::cout << "\n" << std::endl;
    find_chebyshev_center(_mesh, hfs, TetLoop::CHEBY_THRESHOLD, new_pos, true);

    std::cout << "Comparison between star and mesh:\n-Star\n\t-Center = " <<
                 _star.center_ << "\n\t-Bounds = "<<
                 iterableToString<std::set<HalfFaceHandle>>(_star.bounds_) <<
                 " (total= "<< _star.bounds_.size()<< ")"
              << "\n\t-Tets = " <<
                 iterableToString<std::set<CellHandle>>(_star.tets_) << " (total= "<< _star.tets_.size() << ")"
              << "\n-Mesh\n\t-Center = "<<
                 new_pos << "\n\t-Bounds = "<<
                 iterableToString<std::set<HalfFaceHandle>>(hfs) << " (total= "<< hfs.size()<< ")"
              << "\n\t-Tets = " <<
                 iterableToString<std::pair<CellIter, CellIter>>(_mesh.cells())
                                                             << " (total= "<< _mesh.n_logical_cells() << ")"
                                                             << std::endl;
}



bool TetLoop::checkCenter(Star _star){
    bool pass(false);
    double treshold = 2;
    double dist = 0;
    double averageDist = 100;
    for(auto hfh: _star.bounds_){
        auto faceVertices = mesh_.get_halfface_vertices(hfh);
        auto centroid = (mesh_.vertex(faceVertices[0]) +
                mesh_.vertex(faceVertices[1]) +
                mesh_.vertex(faceVertices[2]))/3.f;
        dist += (centroid - _star.center_).norm();
    }
    averageDist = dist / _star.bounds_.size();
    //    std::cout << "Average dist: " << averageDist << std::endl;
    pass = averageDist < treshold;

    return pass;

}

bool TetLoop::checkStarConditions(Star& _star, CellHandle _lastAdded){
    bool isValid(false);
    bool printDebug(false);
    findCavityBoundary(_star);
    isValid = find_chebyshev_center(mesh_, _star.bounds_, CHEBY_THRESHOLD, _star.center_);
    //    std::cout << "Center: "<< _star.center_ <<"\nBounds: "<< _star.bounds_.size() << std::endl;
    if(useWorldSpace_){
        Star starCopy_wm = _star;
        printIterable(starCopy_wm.bounds_);
        findCavityBoundary(starCopy_wm, true);
        isValid = find_chebyshev_center(world_mesh_, starCopy_wm.bounds_, CHEBY_THRESHOLD, starCopy_wm.center_);
        std::cout << "Valid after world mesh cheby check: "<< isValid << std::endl;
    }
    if(!isValid){
        if(printDebug){
            std::cout << "Cheby check failed, revert to previous values" << std::endl;
        }
        _star.tets_.erase(_lastAdded);
        findCavityBoundary(_star);
        find_chebyshev_center(mesh_, _star.bounds_, CHEBY_THRESHOLD, _star.center_);
    }

    return isValid;

}

void TetLoop::findCavityBoundary(std::set<HalfFaceHandle>& _cavityBoundary){
    _cavityBoundary.clear();
    for(auto hfh: mesh_.halffaces()){
        auto opp = mesh_.opposite_halfface_handle(hfh);
        if(cavityEdge_[hfh] && !cavityEdge_[opp] && opp.is_valid()){
            _cavityBoundary.insert(opp);
        }
        if(!opp.is_valid()){
            std::cout << "Opp non valid ! "<< hfh << std::endl;
        }

    }

}

void TetLoop::findCavityBoundary(Star& _constraint, bool _isWorldMesh){
    _constraint.bounds_.clear();
    for(auto tet: _constraint.tets_){
        if(_isWorldMesh){
            for(auto hfh: world_mesh_.cell_halffaces(tet)){
                auto opp = world_mesh_.opposite_halfface_handle(hfh);
                // cavityEdge_wm[hfh] should always be true, since it is part of a star
                if(cavityEdge_wm_[hfh] && !cavityEdge_wm_[opp] && opp.is_valid()){
                    _constraint.bounds_.insert(opp);
                }
                if(!opp.is_valid()){
                    std::cout << "Opp non valid ! "<< hfh << std::endl;
                }
            }
        }else{
            for(auto hfh: mesh_.cell_halffaces(tet)){
                auto opp = mesh_.opposite_halfface_handle(hfh);
                // cavityEdge_[hfh] should always be true, since it is part of a star
                if(cavityEdge_[hfh] && !cavityEdge_[opp] && opp.is_valid()){
                    _constraint.bounds_.insert(opp);
                }
                if(mesh_.is_deleted(opp)){
                    std::cout << "Opp non valid ! "<< hfh << std::endl;
                }
            }
        }
    }
}

CellHandle TetLoop::findNextCell(Star& _startStar, std::vector<Star>& _galaxy){
    bool printDebug(false);
    CellHandle nextCell(-1);
    HalfFaceHandle bestCandidate(-1);
    std::map<double,HalfFaceHandle> candidates;

    for(auto hfh: _startStar.bounds_){
        double f_star = -1;
        auto faceVertices = mesh_.get_halfface_vertices(hfh);
        auto normal = mesh_.normal(hfh);
        auto a = (mesh_.vertex(faceVertices[0]) +
                mesh_.vertex(faceVertices[1]) +
                mesh_.vertex(faceVertices[2]))/3.f;
        f_star = normal.dot(_startStar.center_) - normal.dot(a);
        candidates[f_star] = hfh;
    }
    int size = candidates.size();
    int id = 0;
    do {
        if(candidates.empty()){
            break;
        }
        auto candidateId = *candidates.rbegin();
        bestCandidate = candidateId.second;
        candidates.erase(std::prev(candidates.end()));
        nextCell = mesh_.incident_cell(bestCandidate);
        if(printDebug){
            std::cout << "Trying candidate "<< bestCandidate << " Cell: "<< nextCell << std::endl;
        }
        // check for star collisions
        for(int i = 0; i < _galaxy.size(); ++i){
            auto& star = _galaxy[i];
            if(std::find(star.tets_.begin(), star.tets_.end(), nextCell) != star.tets_.end()){
                // cell is already in another star
                if(printDebug){
                    std::cout << "Candidate is already in a star" << std::endl;
                }
                nextCell = mesh_.InvalidCellHandle;
                break;
            }
        }
        ++id;
    } while (!nextCell.is_valid() && id < size);

    if(!bestCandidate.is_valid()){
        if(printDebug)
            std::cout << "Best candidate halfface is not valid" << std::endl;
        return nextCell;
    }
    if(printDebug){
        if(nextCell.is_valid()){
            std::cout << "Size of boundary: "<< _startStar.bounds_.size() <<"\nNext best cell is "<< nextCell << std::endl;
        }else{
            std::cout << "No valid cell found" << std::endl;
        }

    }

    return nextCell;
}

void TetLoop::triangle_normal_and_centroid(const TetrahedralMesh& _mesh,
                                           const HalfFaceHandle& hf,
                                           ACG::Vec3d& normal,
                                           ACG::Vec3d& centroid,
                                           bool printDebug){
    auto face_vertices = _mesh.get_halfface_vertices(hf);
    auto c_evi  = _mesh.vertex(face_vertices[1]) - _mesh.vertex(face_vertices[0]);
    auto c_evi1 = _mesh.vertex(face_vertices[2]) - _mesh.vertex(face_vertices[0]);


    normal = (c_evi.cross(c_evi1));

    auto normal_norm = normal.norm();
    if(normal_norm){
        normal /= normal_norm;
    }

    centroid = (_mesh.vertex(face_vertices[0]) +
            _mesh.vertex(face_vertices[1]) +
            _mesh.vertex(face_vertices[2]))/3.f;
    if(printDebug){
        std::cout << "--------- checking halfface "<< hf << " (verts: "<<
                     face_vertices[0]<< " "<<
                     face_vertices[1]<< " " <<
                     face_vertices[2]<<
                     ")\n - normal = "<< normal <<
                     "\n - centroid = "<< centroid<< std::endl;
    }

}

bool TetLoop::find_chebyshev_center(const TetrahedralMesh& mesh,
                                    const std::set<HalfFaceHandle>& constraint_hfs,
                                    const double radius_lower_bound,
                                    ACG::Vec3d& new_position,
                                    bool printDebug){

    new_position = {0,0,0};

    using namespace operations_research;
    std::unique_ptr<MPSolver> solver(MPSolver::CreateSolver("GLOP"));

    const double infinity = solver->infinity();


    MPVariable* const x = solver->MakeNumVar(-infinity, infinity, "x");
    MPVariable* const y = solver->MakeNumVar(-infinity, infinity, "y");
    MPVariable* const z = solver->MakeNumVar(-infinity, infinity, "z");
    MPVariable* const r = solver->MakeNumVar(radius_lower_bound, infinity, "r");
    //std::cout << "Number of variables = " << solver->NumVariables();

    for(auto hf: constraint_hfs){

        ACG::Vec3d normal, centroid;

        triangle_normal_and_centroid(mesh, hf,
                                     normal, centroid);

        auto b = normal.dot(centroid);
        if(printDebug){
            std::cout<<" ------- checking halfface "<<mesh.get_halfface_vertices(hf)<<std::endl;
            std::cout<<"  -   normal = "<<normal<<std::endl;
            std::cout<<"  - centroid = "<<centroid<<std::endl;
            std::cout<<"  -        b = "<<b<<std::endl;
        }


        //Ai^T . x + |Ai| . r <= Bi

        MPConstraint* const c = solver->MakeRowConstraint(-infinity, b);
        c->SetCoefficient(x, normal[0]);
        c->SetCoefficient(y, normal[1]);
        c->SetCoefficient(z, normal[2]);

        c->SetCoefficient(r, norm(normal));
    }

    if(printDebug){
        std::cout<<" - #or-tools constraints = " << solver->NumConstraints()<<std::endl;
    }

    //r >= 0
    MPObjective* const objective = solver->MutableObjective();
    objective->SetCoefficient(x, 0);
    objective->SetCoefficient(y, 0);
    objective->SetCoefficient(z, 0);
    objective->SetCoefficient(r, 1);
    objective->SetMaximization();


    //solve LP
    double radius(-1);

    new_position = {0,0,0};
    //compute solution
    const MPSolver::ResultStatus result_status = solver->Solve();
    // Check that the problem has an optimal solution.
    if (result_status != MPSolver::OPTIMAL) {
        //if(result_status == MPSolver::NOT_SOLVED){
        if(printDebug){
            std::cout << "Couldn't solve problem"<<std::endl;
        }
        return false;
    }else{
        new_position[0] = x->solution_value();
        new_position[1] = y->solution_value();
        new_position[2] = z->solution_value();

        radius = r->solution_value();
    }
    if(printDebug){

        std::cout<<" --> or-tools solution: "<<new_position<<std::endl;
        std::cout<<" - radius = "<<radius<<std::endl;
    }

    return true;

}

// ----------- ** Smoothing pass ** -----------------------

void TetLoop::smoothing_pass(PriorityQueue& _A, int _iterations){
    std::set<CellHandle> cellsToAdd;
    for(int i = 0; i < _iterations; ++i){
        // V the vertices of the tets in A
        std::vector<VertexHandle> V;
        while(!_A.empty()) {
            Tet top = _A.top();
            _A.pop();
            auto ch = top.cell_handle_;
            for(auto cvh: mesh_.cell_vertices(ch)){
                V.emplace_back(cvh);
            }
        }
        for(auto v: V){
            Smoothing::smooth(mesh_, v);
            for(auto vch: mesh_.vertex_cells(v)){
                if(mesh_.is_deleted(vch)) continue;
                cellsToAdd.insert(vch);
            }
        }
    }
    for(auto ch: cellsToAdd){
        if(mesh_.is_deleted(ch)) continue;
        double quality = QualityEvaluation::evaluate(ch, mesh_);
        _A.push(Tet(ch,quality));
    }
    cleanQualityQueue(_A, mesh_);
}

// ----------- ** Debugging and tools ** -------------------

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
    mesh.garbage_collection();
    OpenMesh::IO::write_mesh(mesh, LOGS_MESH + "mesh_dump" + std::to_string(id++) + ".obj");
}

void TetLoop::cavityMesh3D(Star& _star){
    static int fileId = 0;
    bool printDebug(false);
    std::set<CellHandle> toAdd;
    TetrahedralMesh copy;
    std::set<VertexHandle> vertices;
    std::map<VertexHandle,VertexHandle> table;
    for(auto ch: _star.tets_){
        toAdd.insert(ch);
    }
    for(auto ch: toAdd){
        auto verts = mesh_.get_cell_vertices(ch);
        for(auto vh: verts){
            if(vertices.find(vh) == vertices.end()){
                table[vh] = copy.add_vertex(mesh_.vertex(vh));
                vertices.insert(vh);
            }
        }
    }
    for(auto ch: toAdd){
        auto verts = mesh_.get_cell_vertices(ch);
        auto c = copy.add_cell(
                    table[verts[0]],
                table[verts[1]],
                table[verts[2]],
                table[verts[3]]);
        std::cout<<" added cell "<<c<<": "<<iterableToString<std::vector<VertexHandle>>(copy.get_cell_vertices(c))<<std::endl;
    }
    if(printDebug){
        std::cout << "Copy with:\n\t-"<<
                     copy.n_cells() << " cells,\n\t-"<<
                     copy.n_faces() << " faces,\n\t-"<<
                     copy.n_edges() << " edges,\n\t-"<<
                     copy.n_vertices()<< " vertices"
                  << std::endl;
    }
    copy.collect_garbage();

    std::cout<<" copy #vertices: "<<copy.n_vertices()<<std::endl;
    std::cout<<" copy cells: "<<std::endl;
    for(auto c: copy.cells()){
        std::cout<<" - "<<c<<": "<<iterableToString<std::vector<VertexHandle>>(copy.get_cell_vertices(c))<<std::endl;
    }

    OpenVolumeMesh::IO::FileManager fm;
    std::string name = "mesh_dump" + std::to_string(fileId++) + "_3D.ovm";
    fm.writeFile(LOGS_MESH + name, copy);
    std::cout << "Created file "<< name << std::endl;
    //    compareStarWithMesh(_star, copy);
    recoverBadStar(_star, name);
}

/**
 * @brief Prints any iterable element with cout
 * @param _toPrint: Iterable to print
 * @param _eol: Should the print end with a line break (optional, default to false)
 */
template <typename T>void TetLoop::printIterable(T _toPrint, bool _eol){
    for(auto elem: _toPrint){
        std::cout << elem << ", ";
    }

    if(_eol){
        std::cout << std::endl;
    }
}

/**
 * @brief Converts an iterable to a string
 * @param _toPrint: Iterable to convert
 * @return A string of form: "x_0, x_1, ..., x_n"
 */
template <typename T> std::string TetLoop::iterableToString(T _toPrint){
    std::string result("");
    for(auto elem: _toPrint){
        result.append(std::to_string(elem.idx()) + ", ");
    }
    result = result.substr(0, result.size()-2);
    return result;
}

/**
 * @brief Converts a vector to a string
 * @param _toPrint: Iterable to convert
 * @return A string of form: "x_0, x_1, ..., x_n"
 */
template <typename T> std::string TetLoop::vectorToString(std::vector<T> _toPrint){
    std::string result("");
    for(auto elem: _toPrint){
        result.append(std::to_string(elem.idx()) + ", ");
    }
    result = result.substr(0, result.size()-2);
    return result;
}

/**
 * @brief Converts a set to a string
 * @param _toPrint: Iterable to convert
 * @return A string of form: "x_0, x_1, ..., x_n"
 */
template <typename T> std::string TetLoop::setToString(std::set<T> _toPrint){
    std::string result("");
    for(auto elem: _toPrint){
        result.append(std::to_string(elem.idx()) + ", ");
    }
    result = result.substr(0, result.size()-2);
    return result;
}

void TetLoop::cleanMesh(TetrahedralMesh& _mesh, bool _keepBoundary){
    bool printDebug(false);
    std::vector<FaceHandle> f_to_del;
    for(auto fh: _mesh.faces()){
        int adjacentCells = 0;
        for(auto f_ch:_mesh.face_cells(fh)){
            ++adjacentCells;
        }
        if(adjacentCells == 0){
            if(printDebug){
                std::cout << "Low valence !!! Face: "<< fh << std::endl;
            }
            f_to_del.push_back(fh);
        }
    }
    for(auto fh: f_to_del){
        _mesh.delete_face(fh);
    }

    std::vector<EdgeHandle> e_to_del;
    for(auto eh: _mesh.edges()){
        int adjacentCells = 0;
        for(auto f_ch:_mesh.edge_cells(eh)){
            ++adjacentCells;
        }
        if(adjacentCells == 0){
            if(printDebug){
                std::cout << "Low valence !!! Edge: "<< eh << std::endl;
            }
            e_to_del.push_back(eh);
        }
    }
    for(auto eh: e_to_del){
        _mesh.delete_edge(eh);
    }

    std::vector<VertexHandle> v_to_del;
    for(auto vh: _mesh.vertices()){
        int adjacentCells = 0;
        for(auto f_ch:_mesh.vertex_cells(vh)){
            ++adjacentCells;
        }
        if(adjacentCells == 0){
            if(_keepBoundary && _mesh.is_boundary(vh)) continue;
            if(printDebug){
                std::cout << "Low valence !!! Vertex: "<< vh << std::endl;
            }
            v_to_del.push_back(vh);
        }
    }
    for(auto vh: v_to_del){
        _mesh.delete_vertex(vh);
    }
}

template <typename T>
void TetLoop::clearDuplicates(std::vector<T>& _vector){
    bool printDebug(false);
    std::set<T> set;
    std::vector<T> noDupl;
    for(auto e: _vector){
        auto inserted = set.insert(e);
        if(inserted.second){
            noDupl.push_back(e);
        }
    }
    if(printDebug){
        std::cout << "Size change with clean: "<< _vector.size() << " to "<< noDupl.size() << std::endl;
    }
    _vector = noDupl;
}

void TetLoop::cleanQualityQueue(TetLoop::PriorityQueue& _queue, TetrahedralMesh& _mesh){
    PriorityQueue copy;
    std::set<CellHandle> unique_handles;
    int size = _queue.size();
    bool printDebug(false);
    while(!_queue.empty()){
        auto top = _queue.top();
        _queue.pop();
        auto inserted = unique_handles.insert(top.cell_handle_);
        if(inserted.second && !_mesh.is_deleted(top.cell_handle_)){
            copy.push(top);
        }
    }
    if(printDebug){
        std::cout << "Queue size change with clean: "
                  << size << " to "<< copy.size() << std::endl;
    }
    _queue = copy;
}

void TetLoop::displayIterationTime(std::chrono::system_clock::time_point& _begin, std::string _name){
    auto end = std::chrono::high_resolution_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - _begin);
    long count = elapsed.count();
    int minutes,seconds,ms;
    // less than a second
    if(count < 1000){
        ms = count;
        std::cout << _name <<" done in " << ms << " ms" <<std::endl;
        return;
    }
    // less than a minute
    if(count < 60000){
        seconds = count/1000;
        ms = count % 1000;
        std::cout << _name <<" done in " << seconds << " s " << ms << " ms" <<std::endl;
        return;
    }
    minutes = count / 60000;
    seconds = count % 60000 / 1000;
    ms = count % 60000 % 1000;
    std::cout << _name <<" done in " << minutes << " min " << seconds << " s " << ms << " ms"<<std::endl;
}

bool TetLoop::validateWorldMesh(TetrahedralMesh& _world_mesh){
    bool valid(false);
    PriorityQueue queue;
    computeQuality(queue, _world_mesh);
    if(queue.top().quality_ > 10e-4){
        valid = true;
    }

    return valid;

}
double TetLoop::computeQualityDelta(double _before){
    bool printDebug(false);
    double delta(0.);
    computeQuality();
    double after = quality_queue_.top().quality_;
    delta = after - _before;
    if(printDebug){
        std::cout << "Quality before: "<< _before
                  << "\t after: "<< after
                  << "\t delta = "<< delta<< std::endl;
    }
    return delta;
}

void TetLoop::addToStats(Stats::StatType _statName, double _quality){
    switch(_statName){
    case Stats::TIMESTEP:
        stats_.timestep_quality_.push_back(_quality);
        break;
    case Stats::TOPOLOGY:
        stats_.topo_quality_delta_.push_back(_quality);
        break;
    case Stats::CONTRACTION:
        stats_.contra_quality_delta_.push_back(_quality);
        break;
    case Stats::INSERTION:
        stats_.insert_quality_delta_.push_back(_quality);
        break;
    case Stats::SMOOTHING:
        stats_.smooth_quality_delta_.push_back(_quality);
        break;
    }

}

void TetLoop::logStats(Stats _stats, Logger& _logger){
    bool printDebug(false);
    if(printDebug){
        std::cout << "Stats\n\tQuality: ";
        printIterable<std::vector<double>>(_stats.timestep_quality_, true);
        std::cout << "\n\tTopological delta: ";
        printIterable<std::vector<double>>(_stats.topo_quality_delta_, true);
        std::cout << "\n\tContraction delta: ";
        printIterable<std::vector<double>>(_stats.contra_quality_delta_, true);
        std::cout << "\n\tInsertion delta: ";
        printIterable<std::vector<double>>(_stats.insert_quality_delta_, true);
        std::cout << "\n\tSmoothing delta: ";
        printIterable<std::vector<double>>(_stats.smooth_quality_delta_, true);
    }

    size_t iter = 0;
    double quality,topo,contra,insert,smooth;
    for(iter; ; ++iter){
        std::vector<double> line(5);
        if(iter < _stats.timestep_quality_.size()){
            quality = _stats.timestep_quality_[iter];
        }else{
            quality = std::numeric_limits<double>::quiet_NaN();
        }
        if(iter < _stats.topo_quality_delta_.size()){
            topo = _stats.topo_quality_delta_[iter];
        }else{
            topo = std::numeric_limits<double>::quiet_NaN();
        }
        if(iter < _stats.contra_quality_delta_.size()){
            contra = _stats.contra_quality_delta_[iter];
        }else{
            contra = std::numeric_limits<double>::quiet_NaN();
        }
        if(iter < _stats.insert_quality_delta_.size()){
            insert = _stats.insert_quality_delta_[iter];
        }else{
            insert = std::numeric_limits<double>::quiet_NaN();
        }
        if(iter < _stats.smooth_quality_delta_.size()){
            smooth = _stats.smooth_quality_delta_[iter];
        }else{
            smooth = std::numeric_limits<double>::quiet_NaN();
        }
        if(std::isnan(quality) && std::isnan(topo)
                && std::isnan(contra) && std::isnan(insert) && std::isnan(smooth)){
            if(printDebug){
                std::cout << "All elements are NaN, end log" << std::endl;
            }
            _logger.close();
            return;
        }
        line = {quality,topo,contra,insert,smooth};
        _logger.logLine(line,true);
    }

}
