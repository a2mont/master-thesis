#include "TriangleLoop.hh"

void TriangleLoop::loop(int _max_iter){
    for(int i=0; i < _max_iter; ++i){
        PriorityQueue worstTriangles;

        computeQuality();

        Triangle currentWorst = quality_queue_.top();
        while(currentWorst.quality_ < q_min_ && !quality_queue_.empty()){
            worstTriangles.push(currentWorst);
            quality_queue_.pop();
            currentWorst = quality_queue_.top();
        }

        if(worstTriangles.empty()){
            std::cout << "All triangles are of good enough quality\nWorst triangle: " << currentWorst.toString() << std::endl;
        }
        else{
            std::cout << worstTriangles.size() << " triangles of bad quality"
                      << "\nWorst triangle: " << worstTriangles.top().toString()
                      << std::endl;
            if(includeLogs_)
                log(logger_);
            improve_mesh(worstTriangles);
        }
        if(includeLogs_)
            log(logger_, true);
        reset_queue(quality_queue_);
    }
    if(includeLogs_)
        log(timeStepLogger_, true);

}

void TriangleLoop::computeQuality(){
    reset_queue(quality_queue_);
    for(auto fh: mesh_.faces()){
        double quality = QualityEvaluation::evaluate(fh, mesh_,false);
        quality_queue_.push(Triangle(fh, quality));
    }
}
void TriangleLoop::log(Logger* _logger, bool _endOfLine){
    computeQuality();
    _logger->logQuality(quality_queue_.top().quality_);
    if(_endOfLine){
        logger_->nextLine();
    }
}

void TriangleLoop::reset_queue(TriangleLoop::PriorityQueue& _queue){
    PriorityQueue empty;
    std::swap(_queue, empty);
}

void TriangleLoop::improve_mesh(PriorityQueue _badTriangles){
    // B <- set of triangles in M with q < q_min
    // foreach t in B
    // if t still exists and q < q_min improve_tet(M,t,q_min)
    while (!_badTriangles.empty()) {
        auto t = _badTriangles.top();
        _badTriangles.pop();
        if(!mesh_.status(t.face_handle_).deleted()){
            improve_triangle(t);
        }
    }
}

void TriangleLoop::improve_triangle(Triangle _t){
    // A <- t, a triangle from the mesh
    PriorityQueue A;
    A.push(_t);
    bool changed;
    int maxIter;
    for (int i = 0; i < 1; ++i) {
        maxIter = 10;
        changed = false;
        do {
            // A <- topological_pass(A,M)
            changed = topologial_pass(&A);
            if(A.top().quality_ >= q_min_){
                if(includeLogs_)
                    log(logger_,true);
                return;
            }
        } while (changed && maxIter-- > 0);

        if(includeLogs_)
            log(logger_);
        // A <- edge_contraction_pass(A,M)
        edge_contraction_pass(&A);
        if(A.top().quality_ >= q_min_){
            if(includeLogs_)
                log(logger_, true);
            return;
        }
        if(includeLogs_)
            log(logger_);
        // A <- insertion_pass(A,M)
        insertion_pass(&A);
        if(A.top().quality_ >= q_min_){
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

bool TriangleLoop::topologial_pass(PriorityQueue* _A){
    //flips,edge removal, multi-face removal operations
    std::cout << "\033[1;36mTopological pass with _A of size: \033[0m" <<  _A->size() << std::endl;
    bool changed = false;
    // E <- edges of the triangles in A
    std::vector<OpenMesh::SmartEdgeHandle> E;
    std::vector<OpenMesh::SmartFaceHandle> newTriangles;

    auto tempMesh = mesh_;
    auto tempA = *_A;

    while(!tempA.empty()){
        auto top = tempA.top();
        tempA.pop();
        auto fh = top.face_handle_;

        for (auto fe_it = tempMesh.fe_cwiter(fh); fe_it.is_valid(); ++fe_it) {
                E.emplace_back(*fe_it);
        }
        // foreach e in E
        // if e still exists attempt to flip e and smooth vertex created
        for(auto e: E){
            if(tempMesh.status(fh).deleted())
                break;

            if (tempMesh.is_flip_ok(e)) {
                std::cout << "Flipped " << e.h0().from() << " , " << e.h0().to() << std::endl;
                tempMesh.flip(e);
                newTriangles.emplace_back(e.h0().face());
                newTriangles.emplace_back(e.h1().face());
                break;
            }
        }
        E.clear();
    }
    // return surviving set in A and the triangles in M altered by the contractions
    for(auto new_fh: newTriangles){
        if (!mesh_.status(new_fh).deleted()) {
            auto qual = QualityEvaluation::evaluate(new_fh,tempMesh);
            tempA.push(Triangle(new_fh, qual));
        }
    }
    if(tempA.top().quality_ > _A->top().quality_){
        *_A = tempA;
        mesh_ = tempMesh;
        changed = true;
        return changed;
    }

    std::cout << "\033[1;33mTopological pass reverted" << std::endl;

    return changed;
}
void TriangleLoop::edge_contraction_pass(PriorityQueue* _A){
    std::cout << "\033[1;36mEdge contraction pass with _A of size: \033[0m" <<  _A->size() << std::endl;
    // E <- edges of the triangles in A
    std::vector<OpenMesh::SmartEdgeHandle> E;
    std::vector<OpenMesh::SmartFaceHandle> newTriangles;

    while(!_A->empty()){
        auto top = _A->top();
        _A->pop();
        auto fh = top.face_handle_;
        E.clear();

        if(mesh_.status(fh).deleted())
            continue;
        for (auto fe_it = mesh_.fe_cwiter(fh); fe_it.is_valid(); ++fe_it) {
                E.emplace_back(*fe_it);
        }
        // foreach e in E
        // if e still exists attempt to contract e and smooth vertex created
        for(auto e: E){
            auto he = e.h0();
            auto from = he.from();
            auto to = he.to();

            if(!mesh_.is_collapse_ok(he) || (mesh_.is_boundary(from) && mesh_.is_boundary(to)))
                continue;
            // if only one vertex is on the boundary, invert the collapse direction
            if(mesh_.is_boundary(from)){
                he = e.h1();
                from = he.from();
                to = he.to();
            }

            // get the faces altered by contraction
            for(auto he_fh: from.faces()){
                newTriangles.emplace_back(he_fh);
            }

            std::cout << "Constraint vertices: " << from.idx() << " "  << to.idx() <<" collapse" << std::endl;
            if(constraint_vhs_.count(from.idx()) > 0){
                constraint_vhs_.erase(from.idx());
                constraint_vhs_[to.idx()] = to.idx();
            }
            mesh_.collapse(he);
            Smoothing::smooth(mesh_, to);
        }
    }
    // return surviving set in A and the triangles in M altered by the contractions
    for(auto new_fh: newTriangles){
        if (!mesh_.status(new_fh).deleted()) {
            _A->push(Triangle(new_fh, QualityEvaluation::evaluate(new_fh,mesh_)));
        }
    }

}



void TriangleLoop::insertion_pass(PriorityQueue* _A){
    std::cout << "\033[1;36mInsertion pass with _A of size: \033[0m" <<  _A->size()  << std::endl;

    std::vector<OpenMesh::SmartFaceHandle> facesContainingP;
    std::vector<OpenMesh::SmartFaceHandle> facesCreated;
    std::map<OpenMesh::SmartVertexHandle, OpenMesh::SmartVertexHandle> borderPairs;
    bool cavityCreated = false;
    // foreach triangle in _A
    while(!_A->empty()){
        // Reset mesh property before pass
        for(auto face: mesh_.faces()){
            mesh_.property(face_visited_, face) = false;
        }
        for(auto heh: mesh_.halfedges()){
            mesh_.property(cavity_edge_, heh) = false;
        }
        borderPairs.clear();
        cavityCreated = false;

        auto top = _A->top();
        _A->pop();
        auto fh = top.face_handle_;

        if(mesh_.status(fh).deleted())
            continue;

        Point p;
        // find the ideal position of new vertex p
        p = mesh_.calc_centroid(fh);
        // Find all existing triangles whose circumscribing circle contains the new point.
        find_faces_with_p(facesContainingP, fh, p);
        for(auto f: facesContainingP){
            if(mesh_.status(f).deleted())
                continue;
            for(auto fheh: f.halfedges()){
                mesh_.property(cavity_edge_, fheh) = true;
                auto opp = fheh.opp();
                if(mesh_.is_boundary(opp)){
                    borderPairs[opp.from()] = opp.to();
                }
            }
            // Update constraint vertices
            for(auto v: f.vertices()){
                if(constraint_vhs_.count(v.idx()>0))
                    constraint_vhs_.erase(v.idx());
            }
        }
        // Delete these triangles; this creates a convex cavity.
        for(auto f: facesContainingP){
            if(mesh_.status(f).deleted())
                continue;
            mesh_.delete_face(f, false);
            cavityCreated = true;
        }
        // If no cavity is created, no need to add a new vertex
        if(!cavityCreated){
            continue;
        }
        auto newVertex = mesh_.add_vertex(p);
        mesh_.set_color(newVertex, ACG::Vec4f(1,1,1,0));
        int tn = 0;
        for(auto heh: mesh_.halfedges()){
            if(!mesh_.status(heh).deleted() && mesh_.property(cavity_edge_, heh)){
                mesh_.property(cavity_edge_, heh) = false;
                //Join the new point to all the vertices on the boundary of the cavity
                auto nfh = mesh_.add_face(newVertex, heh.from(), heh.to());
                std::cout << "Face "<< nfh.idx() << " created with: " << newVertex.idx() << " " << heh.from().idx()
                          << " " << heh.to().idx() <<" Total: "<<++tn <<std::endl;
                facesCreated.emplace_back(nfh);
            }
        }
        for(auto border_vh: borderPairs){
            if(!mesh_.status(border_vh.first).deleted() && !mesh_.status(border_vh.second).deleted()){
                auto nfh = mesh_.add_face(newVertex, border_vh.second, border_vh.first);
                std::cout << "Border face "<< nfh.idx() << " created with: " << newVertex.idx() << " " << border_vh.first.idx()
                          << " " << border_vh.second.idx() <<"Total: "<<++tn <<std::endl;
                facesCreated.emplace_back(nfh);
            }

        }
        Smoothing::smooth(mesh_, newVertex);
    }
    // return surviving set in A and the triangles in M altered by the contractions
    for(auto new_fh: facesCreated){
        if (!mesh_.status(new_fh).deleted()) {
            _A->push(Triangle(new_fh, QualityEvaluation::evaluate(new_fh,mesh_)));
        }
    }
}

void TriangleLoop::find_faces_with_p(std::vector<OpenMesh::SmartFaceHandle>& _list, OpenMesh::SmartFaceHandle _fh, const Point _p){
    if(mesh_.property(face_visited_, _fh) || mesh_.status(_fh).deleted() || !contains_p(_fh, _p))
        return;
    _list.emplace_back(_fh);
    mesh_.property(face_visited_, _fh) = true;
    for(auto ffh: _fh.faces()){
        find_faces_with_p(_list, ffh, _p);
    }

}
// from https://stackoverflow.com/questions/39984709/how-can-i-check-wether-a-point-is-inside-the-circumcircle-of-3-points
bool TriangleLoop::contains_p(OpenMesh::SmartFaceHandle _fh, const Point _p){
    Point points[3];
    int i = 0;
    for(auto vh: _fh.vertices_ccw()){
        points[i++] = mesh_.point(vh);
    }
    Eigen::Matrix3d matrix;
    matrix <<
            points[0][0] - _p[0], points[0][1] - _p[1], pow(points[0][0] - _p[0],2) + pow(points[0][1] - _p[1],2),
            points[1][0] - _p[0], points[1][1] - _p[1], pow(points[1][0] - _p[0],2) + pow(points[1][1] - _p[1],2),
            points[2][0] - _p[0], points[2][1] - _p[1], pow(points[2][0] - _p[0],2) + pow(points[2][1] - _p[1],2);
    if(matrix.determinant() < 0){
        std::cout << _fh.idx() << " contains p: " << _p << std::endl;
        return true;
    }
    return false;
}

void TriangleLoop::smoothing_pass(PriorityQueue* _A, int iterations){
    std::cout << "\033[1;36mSmoothing pass with _A of size: \033[0m" <<  _A->size() << std::endl;
    // V <- the vertices in A
    for(int i = 0; i < iterations; ++i){

        std::vector<OpenMesh::SmartVertexHandle> V;
        for (unsigned long i = 0; i < _A->size(); ++i) {
            Triangle top = _A->top();
            _A->pop();
            auto fh = top.face_handle_;
            for(auto fv_it = mesh_.fv_iter(fh); fv_it.is_valid(); ++fv_it){
                V.emplace_back(*fv_it);
            }
        }

        for(auto v: V){
            Smoothing::smooth(mesh_, v);
            // A <- A U the triangles adjoining v in M
            for(auto vf_it = mesh_.vf_cwiter(v); vf_it.is_valid(); ++vf_it){
                double quality = QualityEvaluation::evaluate(*vf_it, mesh_);
                _A->push(Triangle(*vf_it,quality));
            }
        }
    }
}
