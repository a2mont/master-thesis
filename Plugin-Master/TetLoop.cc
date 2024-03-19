#include "TetLoop.hh"

void TetLoop::loop(int _max_iter){
    auto begin = std::chrono::high_resolution_clock::now();
    double qualityBefore = -std::numeric_limits<double>::infinity();
    static bool collapse(false);
    for(int i=0; i < _max_iter; ++i){
        PriorityQueue worstTets;

        computeQuality();
        if(includeLogs_){
            qualityBefore = quality_queue_.top().quality_;
            double avgQuality(0);
            PriorityQueue q_copy = quality_queue_;
            // register the average quality to stats
            computeQuality();
            while(!q_copy.empty()){
                Tet top = q_copy.top();
                q_copy.pop();
                avgQuality += top.quality_;
            }
            avgQuality /= mesh_.n_logical_cells();

            addToStats(Stats::BEFORE, qualityBefore);
            addToStats(Stats::BEFORE_AVG, avgQuality);
        }


        Tet currentWorst = quality_queue_.top();
        while(currentWorst.quality_ < q_min_ && !quality_queue_.empty()){
            worstTets.push(currentWorst);
            quality_queue_.pop();
            currentWorst = quality_queue_.top();
        }

        if(worstTets.empty()){
            std::cout << "\033[1;32mAll tets are of good enough quality\nWorst tet: \033[0m"
                      << currentWorst.toString() << std::endl;
        }
        else{
            std::cout << worstTets.size() << " tets of bad quality"
                      << "\nWorst tet: " << worstTets.top().toString()
                      << std::endl;
            // subdivision pass on boundary cells
            double delta(0);
//            subdivision_pass(worstTets,delta);
            if(false){
                collapse_surface_pass(worstTets,delta);
                collapse = true;
            }
            improve_mesh(worstTets);
        }
    }
    if(includeLogs_){
        double qualityAfter, avgQuality(0);
        computeQuality();
        qualityAfter = quality_queue_.top().quality_;
        while(!quality_queue_.empty()){
            Tet top = quality_queue_.top();
            quality_queue_.pop();
            avgQuality += top.quality_;
        }
        avgQuality /= mesh_.n_logical_cells();

        addToStats(Stats::AFTER, qualityAfter);
        addToStats(Stats::AFTER_AVG, avgQuality);
    }
    displayIterationTime(begin, "Complete remeshing");
    if(includeLogs_){
        // write to file and reset stats for next timestep
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
    bool topo_active(true), contraction_active(true), insertion_active(true), smoothing_active(true);
    bool changed;
    bool printDebug(false);
    if(printDebug){
        std::cout << "Improve tet: "<< _t.cell_handle_ << std::endl;
    }
    double quality_delta = 0.;
    auto begin = std::chrono::high_resolution_clock::now();
    for (int i = 0; i < 1; ++i) {
        // ensure consistency between spaces, critical error if not valid
        if(useWorldSpace_ && mesh_.n_logical_cells() != world_mesh_.n_logical_cells()){
            std::cout << "World and mesh don't match !" << std::endl;
            std::cout << "Mesh: "<< mesh_.n_logical_cells()
                      << " World: "<< world_mesh_.n_logical_cells() << std::endl;
            exit(EXIT_FAILURE);
        }
        int maxIter = 5;
        changed = false;
        // A <- topological_pass(A,M)
        if(topo_active){
            if(printDebug){
                std::cout << "\033[1;36mTopological pass with A of size: \033[0m" << A.size()<< std::endl;
            }
            auto topo_begin = std::chrono::high_resolution_clock::now();
            double rejected = 0., q_old = quality_queue_.top().quality_;
            do {
                changed = topological_pass(A, quality_delta, rejected);
                if(A.top().quality_ >= q_min_ || A.empty()){
                    if(includeLogs_){
                        addToStats(Stats::TOPOLOGY, quality_delta);
                        addToStats(Stats::TOPOLOGY_REJECT, rejected);
                        quality_delta = 0;
                    }
                    if(printDebug){
                        std::cout << "Exited after topological pass with\n"
                                     "\t-Size of A: "<< A.size() << " -Quality: " <<
                                     A.top().quality_<< std::endl;
                        displayIterationTime(topo_begin, "Topological pass");
                    }
                    return;
                }
            } while (changed && maxIter-- > 0);
            if(includeLogs_){
                addToStats(Stats::TOPOLOGY, quality_delta);
                addToStats(Stats::TOPOLOGY_REJECT, rejected);
                quality_delta = 0;
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
            edge_contraction_pass(A, quality_delta);
            if(A.top().quality_ >= q_min_ || A.empty()){
                if(includeLogs_){
                    addToStats(Stats::CONTRACTION, quality_delta);
                    quality_delta = 0;
                }
                if(printDebug){
                    std::cout << "Exited after contraction pass with\n"
                                 "\t-Size of A: "<< A.size() << " -Quality: " <<A.top().quality_<< std::endl;
                    displayIterationTime(contraction_begin, "Contraction pass");
                }
                return;
            }
            if(includeLogs_){
                addToStats(Stats::CONTRACTION, quality_delta);
                quality_delta = 0;
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
            insertion_pass(A, quality_delta);
            if(A.top().quality_ >= q_min_ || A.empty()){
                if(includeLogs_){
                    addToStats(Stats::INSERTION, quality_delta);
                    quality_delta = 0;
                }
                if(printDebug){
                    std::cout << "Exited after insertion pass with\n"
                                 "\t-Size of A: "<< A.size() << " -Quality: " <<A.top().quality_<< std::endl;
                    displayIterationTime(insertion_begin, "Insertion pass");
                }
                return;
            }
            if(includeLogs_){
                addToStats(Stats::INSERTION, quality_delta);
                quality_delta = 0;
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
            smoothing_pass(A, quality_delta);
            if(A.top().quality_ >= q_min_ || A.empty()){
                if(includeLogs_){
                    addToStats(Stats::SMOOTHING, quality_delta);
                    quality_delta = 0;
                }
                if(printDebug){
                    std::cout << "Exited after smoothing pass with\n"
                                 "-A: "<< A.size() << " , Quality: " <<A.top().quality_<< std::endl;
                    displayIterationTime(smoothing_begin, "Smoothing pass");
                }
                return;
            }
            if(includeLogs_){
                addToStats(Stats::SMOOTHING, quality_delta);
                quality_delta = 0;
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

// ----------- ** Subdivision pass ** ---------------------
void TetLoop::subdivision_pass(PriorityQueue& _A, double& _qualityDelta){
    bool printDebug(false);
    PriorityQueue localQueue;
    double q_old(_A.top().quality_), q_delta(0), q_new; // no default value for q_new
    // récupérer les cells sur la frontière et de mauvaise qualité
    // split l'edge le plus long
    while(!_A.empty()){
        Tet top = _A.top();
        _A.pop();
        CellHandle ch = top.cell_handle_;
        EdgeHandle longest_edge(-1);
        double maxLength = 0;

        // tet is not boundary, add it again to queue
        if(!mesh_.is_boundary(ch)){
            if(printDebug){
                std::cout << "Cell: "<< ch << " not on boundary" << std::endl;
            }
            localQueue.push(top);
            continue;
        }
        // find longest edge on cell
        for(EdgeHandle eh: mesh_.cell_edges(ch)){
            double length = mesh_.length(eh);
            if(length > maxLength){
                longest_edge = eh;
                maxLength = length;
            }
        }
        if(printDebug){
            std::cout << "Longest edge for cell "<< ch
                      << ": "<< longest_edge << " (" << maxLength << ")"<< std::endl;
        }
        // split edge and get new cells
        if(longest_edge.is_valid()){
            VertexHandle added = mesh_.split_edge(longest_edge);
            // same operation on world mesh
            if(useWorldSpace_){
                world_mesh_.split_edge(longest_edge);
            }
            for(CellHandle vch: mesh_.vertex_cells(added)){
                double quality = QualityEvaluation::evaluate(vch, mesh_);
                //TODO:discuss minimal volume or quality ? -------------------------------
                if(quality < q_min_){
                    localQueue.push(Tet(vch, quality));
                }
            }
        }

    }

    q_new = localQueue.top().quality_;
    q_delta = q_new - q_old;
    _qualityDelta += q_delta;

    if(true){
        std::cout << localQueue.size() << " elements added" << std::endl;
        std::cout << "Quality delta: "<< q_delta
                  << " ("<< q_old << " -> " << q_new << ")" << std::endl;
    }
    _A = localQueue;


}

void TetLoop::collapse_surface_pass(PriorityQueue& _A, double& _qualityDelta){
    PriorityQueue localQueue;
    TetrahedralMesh localMesh = mesh_;
    double  treshold(0.1),
            q_old(_A.top().quality_);
    std::vector<std::vector<HalfEdgeHandle>> halfedges_vector;
    std::vector<Tet> tetsToKeep;

    while(!_A.empty()){
        bool collapsable(false);
        Tet top = _A.top();
        _A.pop();
        CellHandle ch = top.cell_handle_;
        if(mesh_.is_deleted(ch) || !mesh_.is_boundary(ch)){
            tetsToKeep.push_back(top);
            continue;
        }
        // get outgoing hes
        for(auto vh: mesh_.get_cell_vertices(ch)){
            if(!mesh_.is_boundary(vh)) continue;
            bool changed(false);
            std::vector<HalfEdgeHandle> halfedges;
            for(auto ohe: mesh_.outgoing_halfedges(vh)){
                halfedges.push_back(ohe);
            }
            for(size_t i = 0; i < halfedges.size()-1; ++i){
                if(mesh_.is_deleted(halfedges[i]) || !mesh_.is_boundary(halfedges[i])) continue;
                VertexHandle from_i(mesh_.from_vertex_handle(halfedges[i])),
                        to_i(mesh_.to_vertex_handle(halfedges[i]));
                Point e0 = mesh_.vertex(to_i) - mesh_.vertex(from_i);
                for(size_t j = i; j < halfedges.size(); ++j){
                    if(mesh_.is_deleted(halfedges[j]) || !mesh_.is_boundary(halfedges[j])) continue;
                    VertexHandle from_j(mesh_.from_vertex_handle(halfedges[j])),
                            to_j(mesh_.to_vertex_handle(halfedges[j]));
                    Point e1 = mesh_.vertex(to_j) - mesh_.vertex(from_j);
                    double product = e0.dot(e1);
                    bool opposite(false);
                    if(product > (-1 - treshold) && product < (-1 + treshold)){
                        opposite = true;
                    }
                    if(!opposite) continue;

                    HalfEdgeHandle heToCollapse(-1);
                    if(link_condition(mesh_, halfedges[i])){
                        heToCollapse = halfedges[i];
                    }else if(link_condition(mesh_, halfedges[j])){
                        heToCollapse = halfedges[j];
                    }
                    if(!heToCollapse.is_valid()) continue;
                    VertexHandle remain = localMesh.collapse_edge(heToCollapse);
                    bool improve(true);
//                    for(auto cell: localMesh.cells()){
//                        double quality = QualityEvaluation::evaluate(cell, localMesh);
//                        if(quality == -std::numeric_limits<double>::infinity()){
//                            improve = false;
//                            break;
//                        }
//                    }
                    if(!improve) continue;
                    remain = mesh_.collapse_edge(heToCollapse);
                    for(auto cell: mesh_.vertex_cells(remain)){
                        double quality = QualityEvaluation::evaluate(cell, mesh_);
                        if(quality < q_min_){
                            localQueue.push(Tet(cell,quality));
                        }
                    }
                    changed = true;
                    collapsable = true;
                    break;
                }
                if(changed) break;
            }
            // if a coollpase occured, go to next cell
            if(changed) break;
        }
        // if no collapse occured, push back cell to queue
        if(!collapsable){
            tetsToKeep.push_back(top);
        }
    }
//    bool changed(false);
//    for(auto halfedges: halfedges_vector){
//        for(size_t i = 0; i < halfedges.size()-1; ++i){
//            if(mesh_.is_deleted(halfedges[i])) continue;
//            VertexHandle from_i(mesh_.from_vertex_handle(halfedges[i])),
//                    to_i(mesh_.to_vertex_handle(halfedges[i]));
//            Point e0 = mesh_.vertex(to_i) - mesh_.vertex(from_i);
//            for(size_t j = i; j < halfedges.size(); ++j){
//                if(mesh_.is_deleted(halfedges[j])) continue;
//                VertexHandle from_j(mesh_.from_vertex_handle(halfedges[j])),
//                        to_j(mesh_.to_vertex_handle(halfedges[j]));
//                Point e1 = mesh_.vertex(to_j) - mesh_.vertex(from_j);
//                double product = e0.dot(e1);
//                bool opposite(false);
//                if(product > -1 - treshold && product < -1 + treshold){
//                    opposite = true;
//                }
//                if(!opposite) continue;

//                bool collapsable(false);
//                if(link_condition(mesh_, halfedges[i])){
//                    collapsable = true;
//                }else if(link_condition(mesh_, halfedges[j])){
//                    collapsable = true;
//                }
//                if(!collapsable) continue;
//                VertexHandle remain = mesh_.collapse_edge(halfedges[i]);
//                for(auto cell: mesh_.vertex_cells(remain)){
//                    double quality = QualityEvaluation::evaluate(cell, mesh_);
//                    if(quality < q_min_){
//                        localQueue.push(Tet(cell,quality));
//                    }
//                }
//                changed = true;
//                break;
//            }
//            if(changed) break;
//        }
//    }
    for(auto& tet: tetsToKeep){
        if(mesh_.is_deleted(tet.cell_handle_)) continue;
        localQueue.push(tet);
    }

    cleanQualityQueue(localQueue,mesh_);
    if(localQueue.empty()){
        return;
    }

    double q_new(localQueue.top().quality_), q_delta(q_new - q_old);
    _qualityDelta = q_delta;
    if(true){
        std::cout << localQueue.size() << " elements added" << std::endl;
        std::cout << "Quality delta: "<< q_delta
                  << " ("<< q_old << " -> " << q_new << ")" << std::endl;
    }
    _A = localQueue;
}


// ----------- ** Topological pass ** ---------------------

bool TetLoop::topological_pass(PriorityQueue& _A, double& _qualityDelta, double& _rejectedTotal){
    //flips,edge removal, multi-face removal operations
    bool changed, edge, face;
    changed = edge = face = false;
    std::vector<CellHandle> cellsAdded;
    std::vector<int> counter({0,0,0});
    double sizeA = _A.size();
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
        if(mesh_.is_deleted(ch)) continue;

        std::vector<EdgeHandle> edges;
        std::vector<FaceHandle> faces;
        for(auto c_eh: mesh_.cell_edges(ch)){
            edges.push_back(c_eh);
        }
        auto edge_begin = std::chrono::high_resolution_clock::now();
        for(auto e: edges){
            if(!mesh_.is_deleted(e)){
                edge = edgeRemoval(e, cellsAdded, _qualityDelta, _rejectedTotal, false);
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
                face = faceRemoval(f, cellsAdded, counter, _qualityDelta, _rejectedTotal, false);
                changed = changed || face;
                total_face = face ? ++total_face : total_face;
            }
        }

//        displayIterationTime(face_begin, "Face part");
        // if nothing happened, we push the cell to the queue again
        if(!changed && !mesh_.is_deleted(ch) && ch.is_valid()){
            cellsAdded.push_back(ch);
        }
    }
    if(total_edge == 0 && total_face == 0 && cellsAdded.size() > sizeA){
        std::cout << "Incorrect amount added" << std::endl;
        exit(EXIT_FAILURE);
    }
//    if(total_edge > 0){
//        std::cout << "Total changes edges: " << total_edge << std::endl;
//    }
//    if(total_face > 0){
//        std::cout << "Method used: \n\t2-3: "<< counter[0]
//                  << "\tMulti: "<< counter[1]<< "\n out of "<< total_face << " changes "
//                  << std::endl;
//    }


    clearDuplicates<CellHandle>(cellsAdded);
    for(auto ch: cellsAdded){
        if(mesh_.is_deleted(ch) || !ch.is_valid()){
            continue;
        }
        double quality = QualityEvaluation::evaluate(ch, mesh_);
        Tet toAdd(ch, quality);
        _A.push(toAdd);
    }
    return changed;
}

// Call when we don't track the cells added by the process (e.g. Tests.cc)
bool TetLoop::edgeRemoval(EdgeHandle _eh,
                          double& _qualityDelta,
                          double& _rejectedTotal,
                          bool _verbose){
    std::vector<CellHandle> ignore;
    return edgeRemoval(_eh, ignore, _qualityDelta, _rejectedTotal, _verbose);
}

bool TetLoop::edgeRemoval(EdgeHandle _eh,
                          std::vector<CellHandle>& _cellsAdded,
                          double& _qualityDelta,
                          double& _rejectedTotal,
                          bool _verbose){
    bool changed = false;
    if(mesh_.is_deleted(_eh) || mesh_.is_boundary(_eh)){
        if(_verbose){
            std::cout << "Edge removal on a deleted/boundary edge" << std::endl;
        }
        return changed;
    }
    TetrahedralMesh wm_backup;
    if(useWorldSpace_){
        wm_backup = world_mesh_;
    }

    // let I the set of tets including eh
    std::set<CellHandle> tetsIncludingEdge;
    // let R the edge ring around eh
    std::set<VertexHandle> oneRingVertices;
    std::vector<HalfEdgeHandle> newHalfEdges;
    PriorityQueue localQueue;
    for(auto ch: mesh_.edge_cells(_eh)){
        if(!ch.is_valid()) continue;
        localQueue.push(Tet(ch,QualityEvaluation::evaluate(ch,mesh_)));
    }

    HalfEdgeHandle toCollapse(-1);
    VertexHandle targetCollapse(-1);

    auto from = mesh_.edge(_eh).from_vertex();
    auto to = mesh_.edge(_eh).to_vertex();

    double q_old = localQueue.top().quality_;
    double maxQuality = q_old;

    for(auto cell: mesh_.edge_cells(_eh)){
        tetsIncludingEdge.insert(cell);
        for(auto vh: mesh_.cell_vertices(cell)){
            if(vh.idx() != from.idx() && vh.idx() != to.idx()){
                oneRingVertices.insert(vh);
            }
        }
    }

    // Empirical results from paper, edge removal for m > 7 rarely improves quality
    if(oneRingVertices.size() > 7){
        return changed = false;
    }

    TetrahedralMesh sub_mesh;

    std::set<VertexHandle> vertices;
    // contains ids of vertice in form table[originalVH] = copyVh
    // !! not valid after mesh changes !!
    std::map<VertexHandle,VertexHandle> meshToCopyVhs;
    // contains ids of vertice in form table[copyVH] = originalVH
    // !! not valid after mesh changes !!
    std::map<VertexHandle,VertexHandle> copyToMeshVhs;
    for(auto ch: mesh_.edge_cells(_eh)){
        if(!ch.is_valid()) continue;
        auto verts = mesh_.get_cell_vertices(ch);
        for(auto vh: verts){
            if(vertices.find(vh) == vertices.end()){
                VertexHandle addedVh = sub_mesh.add_vertex(mesh_.vertex(vh));
                meshToCopyVhs[vh] = addedVh;
                copyToMeshVhs[addedVh] = vh;
                vertices.insert(vh);
            }
        }
    }
    for(auto ch: mesh_.edge_cells(_eh)){
        if(!ch.is_valid()) continue;
        auto verts = mesh_.get_cell_vertices(ch);
        sub_mesh.add_cell(
                meshToCopyVhs[verts[0]],
                meshToCopyVhs[verts[1]],
                meshToCopyVhs[verts[2]],
                meshToCopyVhs[verts[3]]);
    }
    sub_mesh.collect_garbage();
    auto tempEdge = sub_mesh.find_halfedge(meshToCopyVhs[from],meshToCopyVhs[to]);
    VertexHandle newVertex = sub_mesh.split_edge(tempEdge);
    // find the best halfedge to collapse
    for(auto heh: sub_mesh.outgoing_halfedges(newVertex)){
        if(sub_mesh.to_vertex_handle(heh) == meshToCopyVhs[from]
                || sub_mesh.to_vertex_handle(heh) == meshToCopyVhs[to]) continue;
        if(!link_condition(sub_mesh, heh)) continue;
        auto sub_sub_mesh = sub_mesh;
        auto sub_sub_newVh = sub_sub_mesh.collapse_edge(heh);
        std::set<CellHandle> sub_sub_adjacentCells;
        for(auto heh: sub_sub_mesh.incoming_halfedges(sub_sub_newVh)){
            if(oneRingVertices.find(
                        copyToMeshVhs[sub_sub_mesh.from_vertex_handle(heh)]) == oneRingVertices.end()){
                continue;
            }
            for(auto he_ch: sub_sub_mesh.halfedge_cells(heh)){
                sub_sub_adjacentCells.insert(he_ch);
            }
        }
        // Evaluate the quality of the remaining neighbouring cells to the edge
        computeQuality(localQueue,sub_sub_mesh, sub_sub_adjacentCells);

        if(!localQueue.empty() && localQueue.top().quality_ > maxQuality){
            maxQuality = localQueue.top().quality_;
            toCollapse = heh;
            targetCollapse = sub_mesh.to_vertex_handle(toCollapse);
        }
    }

    // collapse the best halfedge on the real mesh
    if(maxQuality > q_old){
        if(!toCollapse.is_valid() || !targetCollapse.is_valid()){
            std::cout << "Invalid !" << std::endl;
            exit(EXIT_FAILURE);
        }
        if(useWorldSpace_){
            TetrahedralMesh wm_copy = world_mesh_;
            VertexHandle added_wm = world_mesh_.split_edge(_eh);
            HalfEdgeHandle sub_to_wm_he =
                    world_mesh_.find_halfedge(added_wm, copyToMeshVhs[targetCollapse]);
            if(!link_condition(world_mesh_, sub_to_wm_he)){
                if(_verbose){
                    std::cout << "Edge removal aborted in world space!!! (cause: link condition)"
                              << std::endl;
                }
                //reset world mesh
                world_mesh_ = wm_copy;
                // missed opportunity for improvement
                _rejectedTotal += -(maxQuality - q_old);
                return changed = false;
            }
            VertexHandle remain_wm = world_mesh_.collapse_edge(sub_to_wm_he);
            Smoothing::smooth(world_mesh_,remain_wm);
            bool valid = validateWorldMesh(world_mesh_);
            if(!valid){
                if(_verbose){
                    std::cout << "Edge removal aborted in world space!!! (cause: invalid)" << std::endl;
                }
                world_mesh_ = wm_copy;
                // missed opportunity for improvement
                _rejectedTotal += -(maxQuality - q_old);
                return changed = false;
            }
        }
        TetrahedralMesh copy = mesh_;
        auto newVertex = mesh_.split_edge(_eh);
        // translate he to collapse from copy to mesh id
        auto translatedHe = mesh_.find_halfedge(newVertex,
                                                copyToMeshVhs[targetCollapse]);
        if(link_condition(mesh_,translatedHe)){
            VertexHandle newFrom = mesh_.collapse_edge(translatedHe);
            Smoothing::smooth(mesh_, newFrom);
            // find the cells adjacent to the remaining vertex after collapse
            std::set<CellHandle> adjacentCells;
            for(auto heh: mesh_.incoming_halfedges(newFrom)){
                if(oneRingVertices.find(mesh_.from_vertex_handle(heh)) == oneRingVertices.end()){
                    continue;
                }
                for(auto he_ch: mesh_.halfedge_cells(heh)){
                    adjacentCells.insert(he_ch);
                }
            }
            computeQuality(localQueue, mesh_);
            if(localQueue.top().quality_ > q_old){
                for(auto ch: adjacentCells){
                    _cellsAdded.push_back(ch);
                }
                _qualityDelta += maxQuality - q_old;
                changed = true;
                if(_verbose){
                    std::cout << "Tets added: "<< vectorToString(_cellsAdded)
                              << ", changed:"<< changed << std::endl;
                }
            }else{
                if(_verbose){
                    std::cout << "Quality outside submesh did not improve" << std::endl;
                }
                mesh_ = copy;
                if(useWorldSpace_){
                    world_mesh_ = wm_backup;
                }
                changed =  false;
            }
        }else{
            if(_verbose){
                std::cout << "Link condition not valid" << std::endl;
            }
            mesh_ = copy;
            if(useWorldSpace_){
                world_mesh_ = wm_backup;
            }
            changed =  false;
        }

        return changed;
    }else{
        double delta = maxQuality - q_old;
        _rejectedTotal += delta;
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
bool TetLoop::faceRemoval(FaceHandle _fh,
                          double& _qualityDelta,
                          double& _rejectedTotal,
                          bool _verbose){
    std::vector<CellHandle> ignore;
    std::vector<int> ig;
    return faceRemoval(_fh, ignore, ig, _qualityDelta, _rejectedTotal, _verbose);
}

bool TetLoop::faceRemoval(FaceHandle _fh,
                          std::vector<CellHandle>& _cellsAdded,
                          std::vector<int>& _counter,
                          double& _qualityDelta,
                          double& _rejectedTotal,
                          bool _verbose){
    bool changed(false);

    if(!_fh.is_valid() || mesh_.is_deleted(_fh) || mesh_.is_boundary(_fh)){
        if(_verbose){
            std::cout << "Face removal on deleted/invalid/boundary face" << std::endl;
        }
        return changed = false;
    }
    TetrahedralMesh subMesh;
    TetrahedralMesh wm_copy = world_mesh_;
    PriorityQueue localQueue;
    std::set<CellHandle> adj;
    std::set<VertexHandle> vertices;
    // contains ids of vertice in form table[originalVH] = copyVh
    // !! not valid after mesh changes !!
    std::map<VertexHandle,VertexHandle> meshToCopyVhs;
    for(auto ch: mesh_.face_cells(_fh)){
        if(!ch.is_valid()) continue;
        auto verts = mesh_.get_cell_vertices(ch);
        for(auto vh: verts){
            if(vertices.find(vh) == vertices.end()){
                VertexHandle addedVh = subMesh.add_vertex(mesh_.vertex(vh));
                meshToCopyVhs[vh] = addedVh;
                vertices.insert(vh);
            }
        }
    }
    for(auto ch: mesh_.face_cells(_fh)){
        if(!ch.is_valid()) continue;
        auto verts = mesh_.get_cell_vertices(ch);
        subMesh.add_cell(
                meshToCopyVhs[verts[0]],
                meshToCopyVhs[verts[1]],
                meshToCopyVhs[verts[2]],
                meshToCopyVhs[verts[3]]);
    }
    subMesh.collect_garbage();
    computeQuality(localQueue, subMesh);
    /* track the changes of the different remeshings
     * ids:
     *  - 0: 2-3 split
     *  - 1: multiface
     *  - 2: Actual result
    */
    std::vector<std::vector<CellHandle>> addedTets(3);
    int bestId = -1;

    // Quality before any operation
    double q_old = localQueue.top().quality_;
    double q_max = q_old;
    double q_23(-std::numeric_limits<double>::infinity()),
            q_multi(-std::numeric_limits<double>::infinity());
    double delta = 0;
    FaceHandle faceSubMesh(-1);
    for(auto face: subMesh.faces()){
        if(!subMesh.is_boundary(face)){
            faceSubMesh = face;
            break;
        }
    }
    if(faceSubMesh.is_valid()){
        flip23(subMesh,faceSubMesh, addedTets[0]);
    }

    // if no new cells were added, skip
    if(!addedTets[0].empty()){
        // Quality after 2-3 flip
        computeQuality(localQueue, subMesh);
        q_23 = localQueue.top().quality_;
        if(_verbose){
            std::cout << "2-3\nNew: " << q_23 << " Max: "<< q_max << std::endl;
        }
        if(q_23 > q_max){
            delta = q_23 - q_max;
            q_max = q_23;
            changed = true;
            bestId = 0;
            if(_verbose){
                std::cout << "2-3 flip has better quality\n" <<
                             "added tets: "<< addedTets[0].size() << std::endl;
            }
        }
    }
//    tempMesh = mesh_;
//    world_mesh_ = wm_copy;
//    multiFace(tempMesh,_fh, addedTets[1]);
//    // if no new cells were added, skip
//    if(!addedTets[1].empty()){
//        computeQuality<std::vector<CellHandle>>(localQueue, tempMesh, addedTets[1]);
//        q_multi = localQueue.top().quality_;
//        if(_verbose){
//            std::cout << "Multi\nNew: " << q_multi << " Max: "<< q_max << std::endl;
//        }
//        if(q_multi > q_max){
//            delta = q_multi - q_max;
//            q_max = q_multi;
//            changed = true;
//            bestId = 1;
//            if(_verbose){
//                std::cout << "Multi face has better quality\n" <<
//                             "Added tets: "<< addedTets[1].size() << std::endl;
//            }
//        }
//    }
    world_mesh_ = wm_copy;
    if(bestId == -1){
        double maxReached = std::max(q_23,q_multi);
        double delta = maxReached - q_old;
        _rejectedTotal += delta;
        return changed = false;
    }
    ++_counter[bestId];
    _qualityDelta += delta;

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
//    displayIterationTime(begin, "Actual choice");

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
    double q_old = std::min({QualityEvaluation::evaluate(a,faceVertices[0],faceVertices[2],faceVertices[1], _mesh),
                             QualityEvaluation::evaluate(faceVertices[1],faceVertices[0],faceVertices[2],b, _mesh),
                             results_uv.o_,
                             results_vw.o_,
                             results_wu.o_
                            });
    double q_new = std::min({results_uv.n_, results_vw.n_, results_wu.n_});
//    std::cout << "Q_new: "<< q_new << " q_old: "<< q_old <<std::endl;
    if(q_new > q_old){
        flip23(_mesh, _fh, _cellsAdded);
        std::vector<FaceWithChildren> total;
//        std::cout << "Result uv h: "<< results_uv.h_.size() << std::endl;
        total.insert(total.end(), results_uv.h_.begin(), results_uv.h_.end());
//        std::cout << "Result vw h: "<< results_vw.h_.size() << std::endl;
        total.insert(total.end(), results_vw.h_.begin(), results_vw.h_.end());
//        std::cout << "Result wu h: "<< results_wu.h_.size() << std::endl;
        total.insert(total.end(), results_wu.h_.begin(), results_wu.h_.end());
//        std::cout << "Total: "<< total.size() << std::endl;
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
    double q_uw = QualityEvaluation::evaluate(a,b,w,u,_mesh);
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
                for(auto& h: result_uv.h_){
                    g.children_.insert(h.fh_);
                }
            }
            if(!result_vw.h_.empty()){
                for(auto& h: result_vw.h_){
                    g.children_.insert(h.fh_);
                }
            }
            double q_old = std::min({QualityEvaluation::evaluate(a,u,w,v,_mesh),
                                     QualityEvaluation::evaluate(v,u,w,b,_mesh),
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
    for(auto& c: cols){
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
//    std::cout << "Children size: "<< _g.children_.size() << std::endl;
    flip32(_mesh, toRemove, _cellsAdded);
    for(auto h: _g.children_){
        flip32Recurse(_mesh, h, _g.fh_, _cellsAdded);
    }
}

void TetLoop::flip32(TetrahedralMesh& _mesh,
                     EdgeHandle _eh,
                     std::vector<CellHandle>& _cellsAdded){
    bool printDebug(false);
    if(_mesh.is_deleted(_eh)|| _mesh.is_boundary(_eh)){
        if(printDebug){
            std::cout << "3-2 flip on a deleted/boundary edge" << std::endl;
        }
        return;
    }

    auto ends = _mesh.edge_vertices(_eh);
    std::vector<HalfEdgeHandle> hes_toCollapse;
    HalfEdgeHandle toCollapse(-1);
    auto tempBaseMesh = _mesh;
    VertexHandle newVertex(-1);
    // the edge ring around eh
    std::set<VertexHandle> oneRingVertices;
    for(auto cell: mesh_.edge_cells(_eh)){
        for(auto vh: mesh_.cell_vertices(cell)){
            if(vh != ends[0] && vh != ends[1]){
                oneRingVertices.insert(vh);
            }
        }
    }

    if(oneRingVertices.size() != 3){
        if(printDebug){
            std::cout << "Incorrect neighbours for 3-2 flip "<< oneRingVertices.size() << std::endl;
        }
        return;
    }
    newVertex = tempBaseMesh.split_edge(_eh);
    if(!newVertex.is_valid()){
        if(printDebug){
            std::cout << "3-2 flip non valid vertex" << std::endl;
        }
        return;
    }
    for(auto heh: tempBaseMesh.outgoing_halfedges(newVertex)){
        if(tempBaseMesh.to_vertex_handle(heh) != ends[0] &&
                tempBaseMesh.to_vertex_handle(heh) != ends[1]){
            hes_toCollapse.push_back(heh);
        }
    }
    // find the best collapse direction
    PriorityQueue localQueue;
    double maxQuality = -std::numeric_limits<double>::infinity();
    for(auto heh: hes_toCollapse){
        auto tempMesh = tempBaseMesh;
        if(!link_condition(tempMesh, heh)) continue;
        auto tempNew = tempMesh.collapse_edge(heh);
        Smoothing::smooth(tempMesh, tempNew);
        // find the cells adjacent to the remaining vertex after collapse
        std::set<CellHandle> adjacentCells;
        for(auto heh: tempMesh.incoming_halfedges(tempNew)){
            if(oneRingVertices.find(tempMesh.from_vertex_handle(heh)) == oneRingVertices.end()){
                continue;
            }
            for(auto he_ch: tempMesh.halfedge_cells(heh)){
                adjacentCells.insert(he_ch);
            }
        }
        for(auto ch: tempMesh.vertex_cells(tempNew)){
            if(!tempMesh.is_deleted(ch) &&
                    std::find(adjacentCells.begin(), adjacentCells.end(), ch)
                    != adjacentCells.end()){
                computeQuality(localQueue,tempMesh, adjacentCells);
            }
        }
        if(!localQueue.empty() && localQueue.top().quality_ > maxQuality){
            maxQuality = localQueue.top().quality_;
            toCollapse = heh;
        }
        reset_queue(localQueue);
    }

    if(link_condition(tempBaseMesh,toCollapse)){
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
        _mesh.split_edge(_eh);
        auto remain = _mesh.collapse_edge(toCollapse);
        Smoothing::smooth(_mesh,remain);
        // find the cells adjacent to the remaining vertex after collapse
        std::set<CellHandle> adjacentCells;
        for(auto heh: _mesh.incoming_halfedges(remain)){
            if(oneRingVertices.find(_mesh.from_vertex_handle(heh)) == oneRingVertices.end()){
                continue;
            }
            for(auto he_ch: _mesh.halfedge_cells(heh)){
                adjacentCells.insert(he_ch);
            }
        }
        for(auto ch: adjacentCells){
            _cellsAdded.push_back(ch);
        }
    }
}

void TetLoop::flip23(TetrahedralMesh& _mesh,
                     FaceHandle _fh,
                     std::vector<CellHandle>& _cellsAdded){
    bool printDebug(true);
    if(_mesh.is_deleted(_fh) || _mesh.is_boundary(_fh)){
        if(printDebug){
            std::cout << "2-3 flip on a deleted/boundary face" << std::endl;
        }
        return;
    }
    auto hfh = _mesh.face_halffaces(_fh)[0];
    auto toVertex = _mesh.halfface_opposite_vertex(hfh);
    if(!_mesh.is_valid(toVertex)){
        hfh = _mesh.face_halffaces(_fh)[1];
        toVertex = _mesh.halfface_opposite_vertex(hfh);
    }
    if(!_mesh.is_valid(toVertex)){
        if(printDebug){
            std::cout << "No valid vertex found for 2-3 flip" << std::endl;
        }
        return;
    }
    auto hfh_opp = _mesh.opposite_halfface_handle(hfh);
    auto vh_opp = _mesh.halfface_opposite_vertex(hfh_opp);

    auto temp = _mesh;
    auto newVertex = temp.split_face(_fh);
    auto heToCollapse = temp.find_halfedge(newVertex, toVertex);

    if(link_condition(temp, heToCollapse)){
        if(useWorldSpace_){
            auto wm_copy = world_mesh_;
            // update world space
            VertexHandle added_wm = world_mesh_.split_face(_fh);
            world_mesh_.collapse_edge(heToCollapse);
            Smoothing::smooth(world_mesh_, added_wm);
            bool valid = validateWorldMesh(world_mesh_);
            // operation is not valid in world space -> abort
            if(!valid){
                world_mesh_ = wm_copy;
                if(printDebug){
                    std::cout << "2-3 split aborted in world space!!!" << std::endl;
                }
                return;
            }
        }
        VertexHandle added = _mesh.split_face(_fh);
        auto remain = _mesh.collapse_edge(heToCollapse);
        Smoothing::smooth(_mesh, added);
        auto remainingHe = _mesh.find_halfedge(vh_opp, remain);
        for(auto ch: _mesh.halfedge_cells(remainingHe)){
            if(!_mesh.is_deleted(ch)){
                _cellsAdded.push_back(ch);
            }

        }

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

void TetLoop::edge_contraction_pass(PriorityQueue& _A, double& _qualityDelta){
    bool printDebug(false);
    std::vector<EdgeHandle> edges;
    std::vector<CellHandle> newTets;
    PriorityQueue tempA = _A;
    TetrahedralMesh mesh_copy = mesh_;
    double rejected(0);

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
            auto remain = contractEdge(e, newTets, _qualityDelta, rejected);
            changed = changed || remain.is_valid();
        }
        // if nothing happenend to the tet, push back onto queue
        if(!changed){
            newTets.push_back(ch);
        }

    }

    addToStats(Stats::CONTRACTION_REJECT, rejected);

    if(newTets.empty()){
        cleanQualityQueue(tempA, mesh_);
        _A = tempA;
        if(printDebug){
            std::cout << "new cells empty" << std::endl;
        }
        return;
    }

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
        std::cout << "Total: "<< newTets.size()
                  << " vs duplicates/deleted removed: "
                  << newTets_no_dupl.size() << std::endl;
    }
}

VertexHandle TetLoop::contractEdge(EdgeHandle _eh,
                                   std::vector<CellHandle>& _tetsAltered,
                                   double& _qualityDelta,
                                   double& _rejectTotal){
    bool printDebug(false);
    VertexHandle remain = mesh_.InvalidVertexHandle;
    std::vector<CellHandle> cellsAround;
    for(auto e_ch: mesh_.edge_cells(_eh)){
        cellsAround.push_back(e_ch);
    }
    PriorityQueue queue;
    computeQuality<std::vector<CellHandle>>(queue,mesh_,cellsAround);
    double qualityBefore = queue.top().quality_;

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
        std::vector<CellHandle> around;
        for(auto v_ch: temp.vertex_cells(tempRemain)){
            around.push_back(v_ch);
        }
        PriorityQueue tempQueue;
        computeQuality<std::vector<CellHandle>>(tempQueue, temp, around);
        // accept only operations that improve quality
        if(tempQueue.top().quality_ > qualityBefore){
            if(useWorldSpace_){
                auto wm_copy = world_mesh_;
                wm_copy.collapse_edge(he);
                bool valid = validateWorldMesh(wm_copy);
                // operation is not valid in world space -> abort
                if(!valid){
                    if(printDebug){
                        std::cout << "Edge contraction aborted in world space!!!" << std::endl;
                    }
                    _rejectTotal += -(tempQueue.top().quality_ - qualityBefore);
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
            _qualityDelta += tempQueue.top().quality_ - qualityBefore;
            // get the tets around new point
            for(auto v_ch: mesh_.vertex_cells(remain)){
                _tetsAltered.push_back(v_ch);
            }
        }else{
            _rejectTotal += tempQueue.top().quality_ - qualityBefore;
        }
    }

    return remain;
}

// ----------- ** Insertion pass ** -----------------------


void TetLoop::insertion_pass(PriorityQueue& _A, double& _qualityDelta){
    bool saveMesh(false);
    bool printDebug(false);

    std::vector<CellHandle> cellsToAdd;
    std::vector<Star> galaxy;
    std::vector<Star> galaxy_wm;
    CellHandle lastAdded(-1);

    double delta(0), reject(0);
    const int maxStarSize = 5;

    auto copy = mesh_;
    auto wm_copy = world_mesh_;
    auto wm_temp = world_mesh_;
    // copy of A to revert to
    PriorityQueue copyA = _A;
    // used in the topo pass
    PriorityQueue topo_queue;
    // used to calculate the quality delta for each star
    PriorityQueue localQueue;
    std::vector<double> starQualities;
    computeQuality();
    double q_old = quality_queue_.top().quality_;
    //    for each inverted tetrahedron 𝑐 ∈ 𝐶
    while(!_A.empty()){
        for(auto hfh: mesh_.halffaces()){
            cavityEdge_[hfh] = false;
            if(useWorldSpace_){
                cavityEdge_wm_[hfh] = false;
            }
        }
        auto top = _A.top();
        _A.pop();
        auto ch = top.cell_handle_;
        if(mesh_.is_deleted(ch)) continue;
        double starCount(0);
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
        starCount = 1;
        for(auto hfh: mesh_.cell_halffaces(ch)){
            auto opp = mesh_.opposite_halfface_handle(hfh);
            cavityEdge_[hfh] = true;
            if(useWorldSpace_){
                cavityEdge_wm_[hfh] = true;
            }
            if(opp.is_valid()){
                newStar.bounds_.insert(opp);
            }else{
                std::cout << "Non valid halfface added to star bounds" << std::endl;
            }
        }
        if(newStar.bounds_.size() == 4){
            find_chebyshev_center(mesh_, newStar.bounds_, CHEBY_THRESHOLD, newStar.center_);
        }else{
            std::cout <<
                         "Error on initial cavity boundary size !!\n\tExpected: 4  Value: "
                      << newStar.bounds_.size()<< std::endl;
        }
        // used to check if no neighbour can be added to star
        bool optionsLeft(true);
        do {
            //    Choose next tetrahedron 𝑐∗ by heuristic
            auto next = findNextCell(newStar, galaxy, optionsLeft);
            if(next.is_valid()){
                //    𝑆 ∗ = ⟨𝑆 ∗ ∪ {𝑐 ∗ } ⟩
                ++starCount;
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
        } while (checkStarConditions(newStar, lastAdded) && starCount < maxStarSize && optionsLeft);
        //    G = G ∪ {𝑆 ∗ }
        galaxy.push_back(newStar);
        if(useWorldSpace_){
            Star wmStar = newStar;
            find_chebyshev_center(wm_temp, wmStar.bounds_, CHEBY_THRESHOLD, wmStar.center_);
            galaxy_wm.push_back(wmStar);
        }
    }
    if(useWorldSpace_){
        if(galaxy.size() != galaxy_wm.size()){
            std::cout << "Galaxy size missmatch !!" << std::endl;
            std::cout << "Galaxy sizes, mesh: "<< galaxy.size()
                      << " world: "<< galaxy_wm.size() << std::endl;
            exit(EXIT_FAILURE);
        }
    }
    for(size_t i = 0; i < galaxy.size();++i){
        auto& star = galaxy[i];
        if(!checkCenter(star)){
            if (printDebug) {
                bool cheb_test = find_chebyshev_center(mesh_, star.bounds_, CHEBY_THRESHOLD, star.center_);
                std::cout << "Center of star check failed with cheb test: "<< cheb_test << std::endl;
            }
            if(saveMesh){
                cavityMesh3D(star);
            }
            // If center fails (e.g. single cell is too small to find a cheby center), cancel
            continue;
        }
        for(auto hfh: star.bounds_){
            if(mesh_.is_deleted(hfh)){
                std::cout<<" star bounds hf is deleted"<<std::endl;
            }
            // add vertices for cell reconstruction
            auto vertices = mesh_.get_halfface_vertices(hfh);
            star.reconstructionVectors_.push_back(vertices);

        }
        if(printDebug && star.bounds_.size() != star.reconstructionVectors_.size()){
            std::cout << "Not enough vertices triplet to reconstruct star!\n-Bounds size:"
                      << star.bounds_.size()<<
                      " Reconstruct size: "<< star.reconstructionVectors_.size()<< std::endl;
        }

        // save quality of star before operation
        computeQuality<std::set<CellHandle>>(localQueue,mesh_,star.tets_);
        starQualities.push_back(localQueue.top().quality_);

        for(auto ch: star.tets_){
            if(mesh_.is_deleted(ch)){
                std::cout << "Already deleted cell: " << ch << std::endl;
                continue;
            }
            mesh_.delete_cell(ch);
        }
        if(useWorldSpace_){
            Star& wm_star = galaxy_wm[i];
            if(star.tets_.size() != wm_star.tets_.size()){
                std::cout << "Tets size missmatch !!!" << std::endl;
                std::cout << "mesh "<< star.tets_.size() << std::endl;
                std::cout << "text "<< wm_star.tets_.size() << std::endl;
            }
            for(auto ch: wm_star.tets_){
                if(wm_temp.is_deleted(ch)){
                    std::cout << "Already deleted cell (world mesh): " << ch << std::endl;
                    continue;
                }
                wm_temp.delete_cell(ch);
            }
        }
    }
    if(useWorldSpace_){
        cleanMesh(world_mesh_,true);
    }
    // Remplir la galaxy
    // for each star
    // 1. compute center of cheby and add point
    // 2. fill with vertices from boundary (reconstruction vertices + newPoint)
    for(size_t i = 0; i < galaxy.size(); ++i){
        std::set<CellHandle> cavityFill;
        auto& star = galaxy[i];
        if(!checkCenter(star)){
            // If center fails (e.g. single cell is too small to find a cheby center), cancel
            continue;
        }
        // copy to revert to in case world space aborts operation
        auto meshCopy = mesh_;
        auto wmCopy = wm_temp;
        VertexHandle center_wm;
        if(useWorldSpace_){
            Star& star_wm = galaxy_wm[i];
            center_wm = wm_temp.add_vertex(star_wm.center_);
        }
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

            Smoothing::smooth(mesh_, newVertex);
            if(useWorldSpace_){
                wm_temp.add_cell(vertices[0], vertices[2], vertices[1], center_wm);
            }
            //            std::cout << "Added cell "<< added.idx() << " valid = "<< added.is_valid() << std::endl;
            if(added.is_valid()){
                cavityFill.insert(added);
            }
        }
        // checks if the correct number of tets has been added
        size_t added_size = cavityFill.size();
        if(added_size != star.reconstructionVectors_.size()){
            std::cout << "Added size: "<< added_size <<
                         " Recons vector size: "<< star.reconstructionVectors_.size() << std::endl;
            exit(EXIT_FAILURE);
        }
        if(useWorldSpace_){
            world_mesh_ = wm_temp;
        }
        // topological pass on the new cells
        int attempts = 1;
        bool changed(false);
        double ignoreDelta(0);
        for(auto ch: cavityFill){
            double quality = QualityEvaluation::evaluate(ch,mesh_);
            topo_queue.push(Tet(ch,quality));
        }
        do{
            if(printDebug){
                std::cout << "Attempt: "<< attempts << std::endl;
                std::cout << "Worst element before insertion topo pass: "<< topo_queue.top().quality_
                          << std::endl;
            }
            changed = topological_pass(topo_queue, ignoreDelta, ignoreDelta);
            if(printDebug){
                std::cout << "Worst element after insertion topo pass: "<< topo_queue.top().quality_
                          << std::endl;
            }
        }while(--attempts > 0 && changed);
        // copy to add to local queue
        PriorityQueue topo_queue_copy = topo_queue;
        while(!topo_queue_copy.empty()){
            Tet top = topo_queue_copy.top();
            topo_queue_copy.pop();
            localQueue.push(top);
        }
        // get star quality delta
        if(!topo_queue.empty()){
            double qualityAfter = topo_queue.top().quality_;
            delta+= qualityAfter - starQualities[i];
        }
        if(useWorldSpace_){
            bool valid = validateWorldMesh(world_mesh_);
            if(!valid){
                if(delta > 0){
                    reject -= delta;
                }else{
                    reject += delta;
                }
                mesh_ = meshCopy;
                world_mesh_ = wmCopy;
                if(printDebug){
                    std::cout << "Insertion cancelled by world mesh" << std::endl;
                }
                continue;
            }
        }
    }
    // Ensure that no unlinked elements remain
    cleanMesh(mesh_);
    if(useWorldSpace_){
        cleanMesh(world_mesh_);
    }
    // ensure that the mesh quality improves
    if(delta < 1e-6){
        reject = delta;
        addToStats(Stats::INSERTION_REJECT, reject);
        mesh_ = copy;
        cleanQualityQueue(copyA, mesh_);
        _A = copyA;
        if(useWorldSpace_){
            world_mesh_ = wm_copy;
        }
        if(printDebug){
            std::cout << "Quality did not improve, revert (delta = "
                      << delta << ")" << std::endl;
        }
        return;
    }
    computeQuality();
    if(quality_queue_.top().quality_ < q_old){
        delta = quality_queue_.top().quality_ - q_old;
        reject = delta;
        addToStats(Stats::INSERTION_REJECT, reject);
        mesh_ = copy;
        cleanQualityQueue(copyA, mesh_);
        _A = copyA;
        if(useWorldSpace_){
            world_mesh_ = wm_copy;
        }
        if(true){
            std::cout << "Global quality did not improve, revert (delta = "
                      << delta << ")" << std::endl;
        }
        return;
    }

    if(printDebug){
        std::cout << "Improvement "<< delta << std::endl;
    }
    _qualityDelta += delta;
    addToStats(Stats::INSERTION_REJECT, reject);
    _A = localQueue;
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
    find_chebyshev_center(mesh_, _star.bounds_, CHEBY_THRESHOLD, _star.center_);

    std::set<HalfFaceHandle> hfs;
    for(auto hf:_mesh.halffaces()){
        if(_mesh.is_boundary(hf)){
            hfs.insert(hf);
        }
    }
    ACG::Vec3d new_pos;
    find_chebyshev_center(_mesh, hfs, TetLoop::CHEBY_THRESHOLD, new_pos);

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
    if(_star.bounds_.size() != hfs.size()){
        std::cout << "Inequal bound size !!: Star "<< _star.bounds_.size()
                  << " vs Mesh "<< hfs.size() << std::endl;
        std::cout << "Boundary halffaces in star:" << std::endl;
        for(auto hfh: _star.bounds_){
            if (mesh_.is_boundary(hfh)){
                std::cout << hfh <<", ";

            }
        }
        TetrahedralMesh copy = mesh_;

        copy.collect_garbage();

        OpenVolumeMesh::IO::FileManager fm;
        std::string name = "mesh_dump_bounds_3D.ovm";
        fm.writeFile(LOGS_MESH + name, copy);
        std::cout << "Created file "<< name << std::endl;
        exit(EXIT_FAILURE);
    }

}



bool TetLoop::checkCenter(Star _star){
    bool pass(false);
    double treshold = 2.5;
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
//        printIterable(starCopy_wm.bounds_);
        findCavityBoundary(starCopy_wm, true);
        isValid = find_chebyshev_center(world_mesh_, starCopy_wm.bounds_,
                                        CHEBY_THRESHOLD, starCopy_wm.center_);
//        std::cout << "Valid after world mesh cheby check: "<< isValid << std::endl;
    }
    if(!isValid){
        if(printDebug){
            std::cout << "Cheby check failed, revert to previous values" << std::endl;
        }
        for(auto hfh: mesh_.cell_halffaces(_lastAdded)){
            cavityEdge_[hfh] = false;
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

CellHandle TetLoop::findNextCell(Star& _startStar, std::vector<Star>& _galaxy, bool& _resultsLeft){
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
        for(size_t i = 0; i < _galaxy.size(); ++i){
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

    if(!nextCell.is_valid()){
        // we checked every neighbour, none work
        _resultsLeft = false;
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

void TetLoop::smoothing_pass(PriorityQueue& _A, double& _qualityDelta, int _iterations){
    std::set<CellHandle> cellsToAdd;
    PriorityQueue copyA = _A;
    TetrahedralMesh meshCopy = mesh_;
//    computeQuality();
    double reject(0), totalAdded(0);//, q_old(quality_queue_.top().quality_);
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
                double before(0),after(0),delta(0);
                Point posBefore = mesh_.vertex(v);
                for(auto vch: mesh_.vertex_cells(v)){
                    if(cellsToAdd.find(vch) == cellsToAdd.end()){
                        before += QualityEvaluation::evaluate(vch, mesh_);
                    }
                }
                Smoothing::smooth(mesh_, v);
                for(auto vch: mesh_.vertex_cells(v)){
                    if(cellsToAdd.find(vch) == cellsToAdd.end()){
                        after += QualityEvaluation::evaluate(vch, mesh_);
                    }
                }
                delta = after - before;
                //ensure mesh improvement
                if(delta < 1e-6){
                    reject += delta;
                    // reset vertex to pos before smoothing
                    mesh_.set_vertex(v, posBefore);
                    continue;
                }

                _qualityDelta+= delta;
                for(auto vch: mesh_.vertex_cells(v)){
                    cellsToAdd.insert(vch);
                }
            }
        }
        addToStats(Stats::SMOOTHING_REJECT, reject);
        for(auto ch: cellsToAdd){
            if(mesh_.is_deleted(ch)) continue;
            double quality = QualityEvaluation::evaluate(ch, mesh_);
            _A.push(Tet(ch,quality));
        }
        cleanQualityQueue(_A, mesh_);
}



// ----------- ** Debugging and tools ** -------------------

void TetLoop::compute_volumes_and_min_heights(TetrahedralMesh& mesh,
                                     double& min_height,
                                     double& min_vol){

    min_height = std::numeric_limits<double>::max();

    for(auto c: mesh.cells()){
        //std::cout<<" ----------- cell "<<c<<":"<<std::endl;

        auto vol = QualityEvaluation::calculate_volume(mesh,c);

        min_vol = std::min(min_vol, vol);

        //skip degenerate cells
        if(vol == 0){
            continue;
        }

        auto c_vertices = mesh.get_cell_vertices(c);
        auto e1 = (mesh.vertex(c_vertices[1]) - mesh.vertex(c_vertices[0]));
        auto e2 = (mesh.vertex(c_vertices[2]) - mesh.vertex(c_vertices[0]));
        auto e3 = (mesh.vertex(c_vertices[3]) - mesh.vertex(c_vertices[0]));

        auto max_area = 0.5 * std::max((e1.cross(e2)).norm(),
                                       std::max((e1.cross(e3)).norm(),
                                                (e2.cross(e3)).norm()));
        //std::cout<<"  - max area = "<<max_area<<std::endl;

        if(max_area != 0){
            auto cell_min_height = 6.0 * vol / max_area;

            min_height = std::min(min_height, cell_min_height);
        }else{

            min_height = 0;
        }
    }
}

void TetLoop::computeQuality(){
    reset_queue(quality_queue_);
    for(auto c_it = mesh_.cells_begin(); c_it != mesh_.cells_end(); ++c_it){
        if(!c_it->is_valid() || mesh_.is_deleted(*c_it)){
            std::cout <<*c_it <<" is deleted" << std::endl;
            continue;
        }
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

template <typename T>
void TetLoop::computeQuality(TetLoop::PriorityQueue& _queue,
                             TetrahedralMesh& _mesh, const T& _iterableToEvaluate){
    reset_queue(_queue);
    for(auto ch: _iterableToEvaluate){
        // minus to handle the change of quality to dirichlet
        double quality = QualityEvaluation::evaluate(ch, _mesh);
        _queue.push(Tet(ch, quality));
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
    if(printDebug){
        for(auto line: table){
            std::cout << line.first << " -> "<< line.second << std::endl;
        }
    }
    for(auto ch: toAdd){
        auto verts = mesh_.get_cell_vertices(ch);
        auto c = copy.add_cell(
                    table[verts[0]],
                table[verts[1]],
                table[verts[2]],
                table[verts[3]]);
        if(printDebug){
            std::cout<<" added cell "<<c<<": "<<
                       iterableToString<std::vector<VertexHandle>>
                       (copy.get_cell_vertices(c))<<std::endl;
        }
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


    for(auto hfh: _star.bounds_){
        auto vertices = mesh_.get_halfface_vertices(hfh);
        std::vector<VertexHandle> v_vec({table[vertices[0]],table[vertices[1]], table[vertices[2]]});
        HalfFaceHandle hf_copy = copy.find_halfface(v_vec);
        std::cout << "Halfface "<< hfh << " -> "<< hf_copy << std::endl;
    }
    for(auto hfh: copy.halffaces()){
        if(copy.is_boundary(hfh)){
            auto vertices_copy = copy.get_halfface_vertices(hfh);
            std::vector<VertexHandle> vertices;
            HalfFaceHandle key(-1);
            for(auto vh: vertices_copy){
                for(auto pair: table){
                    if(pair.second == vh){
                        vertices.push_back(pair.first);
                    }
                }
            }
            key = mesh_.find_halfface(vertices);
            if(_star.bounds_.find(key) == _star.bounds_.end()){
                std::cout << "Mesh Hface "<< hfh
                          << " is boundary, Star Hface "<< key << " is not on star boundary"<< std::endl;
                std::cout << "Cavity edge property for "<< key <<": "<< cavityEdge_[key] << std::endl;
                auto temp = mesh_;
                auto cell = mesh_.incident_cell(key);
                for(auto ch: temp.cells()){
                    if(ch != cell){
                        temp.delete_cell(ch);
                    }
                }
                temp.collect_garbage();

                OpenVolumeMesh::IO::FileManager fm;
                std::string name = "star_only_boundary_3D.ovm";
                fm.writeFile(LOGS_MESH + name, temp);


            }

        }
    }

    if(printDebug){
        std::cout<<" copy #vertices: "<<copy.n_vertices()<<std::endl;
        std::cout<<" copy cells: "<<std::endl;
        for(auto c: copy.cells()){
            std::cout<<" - "<<c<<": "<<iterableToString<std::vector<VertexHandle>>(copy.get_cell_vertices(c))<<std::endl;
        }
    }

    OpenVolumeMesh::IO::FileManager fm;
    std::string name = "mesh_dump" + std::to_string(fileId++) + "_3D.ovm";
    fm.writeFile(LOGS_MESH + name, copy);
    std::cout << "Created file "<< name << std::endl;
    compareStarWithMesh(_star, copy);
//    recoverBadStar(_star, name);
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

void TetLoop::saveMesh(TetrahedralMesh& _mesh, std::set<CellHandle>& _cellsToKeep, std::string _name){
    TetrahedralMesh copy = _mesh;
    for(auto ch: copy.cells()){
        if(_cellsToKeep.find(ch) == _cellsToKeep.end()){
            copy.delete_cell(ch);
        }
    }
    cleanMesh(copy);
    copy.collect_garbage();
    OpenVolumeMesh::IO::FileManager fm;
    fm.writeFile(LOGS_MESH + _name, copy);
    std::cout << "Created file "<< _name << std::endl;
}

void TetLoop::saveMesh(TetrahedralMesh& _mesh, std::string _name){
    TetrahedralMesh copy = _mesh;
    cleanMesh(copy);
    copy.collect_garbage();
    OpenVolumeMesh::IO::FileManager fm;
    fm.writeFile(LOGS_MESH + _name, copy);
    std::cout << "Created file "<< _name << std::endl;
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

bool TetLoop::validateWorldMesh(TetrahedralMesh& _mesh){
    bool valid(true);
    for(auto ch: _mesh.cells()){
        double quality = QualityEvaluation::evaluate(ch, _mesh);
        if(quality == -std::numeric_limits<double>().infinity()){
            valid = false;
            break;
        }
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
//    if(_quality == -std::numeric_limits<double>::infinity() ||
//            _quality == std::numeric_limits<double>::infinity()){
//        stats_.stat_data_[_statName].push_back(0);
//    } else{
        stats_.stat_data_[_statName].push_back(_quality);
//    }
}

void TetLoop::logStats(Stats& _stats, Logger& _logger){
    bool printDebug(false);
    if(printDebug){
        for(Stats::StatType i = Stats::INIT; i != Stats::LAST; i = Stats::StatType(i+1)){
            Stats::StatType type = static_cast<Stats::StatType>(i);
            std::cout << _stats.types_titles[type] << std::endl;
            printIterable<std::vector<double>>(_stats.stat_data_[type], true);
        }
    }

    size_t iter = 0;
    for(iter; ; ++iter){
        std::vector<double> line(_stats.stat_data_.size());
        for(auto& kv_pair: _stats.stat_data_){
            if(iter < kv_pair.second.size()){
                line[kv_pair.first] = kv_pair.second[iter];
            }else{
                line[kv_pair.first] = std::numeric_limits<double>::quiet_NaN();
            }
        }
        bool complete(false);
        // checks that all stats values are done recording
        for(auto l: line){
            if(!std::isnan(l)){
                complete = false;
                break;
            }
            complete = true;
        }
        if(complete){
            if(printDebug){
                std::cout << "All elements are NaN, end log" << std::endl;
            }
            _logger.close();
            return;
        }
        _logger.logLine(line,true);
    }

}
