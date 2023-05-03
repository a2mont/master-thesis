#include "MainLoop.hh"



void MainLoop::loop(ACG::Vec3d _displacement, int _constraint_vh, bool _verbose,int _max_iter ){

    for(int i=0; i < _max_iter; ++i){
        PriorityQueue worstTriangles;

        vd_.displace(_displacement, _constraint_vh, _verbose);

        for(auto fh: mesh_.faces()){
            double quality = QualityEvaluation::evaluate(fh, mesh_,_verbose);
            quality_queue_.push(Triangle(fh, quality));
        }

        Triangle worstTriangle = quality_queue_.top();
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
        if(!mesh_.status(t.face_handle_).deleted())
            improve_triangle(t);
    }
}

void MainLoop::improve_triangle(Triangle _t){
    // A <- t, a triangle from the mesh
    PriorityQueue A;
    A.push(_t);
    bool changed;
    for (int i = 0; i < 10; ++i) {
        changed = false;
        do {
            // A <- topological_pass(A,M)
//            changed = topologial_pass(&A);
            if(A.top().quality_ >= q_min_)
                return;
        } while (changed);

        // A <- edge_contraction_pass(A,M)
//        edge_contraction_pass(&A);
        if(A.top().quality_ >= q_min_)
            return;

        // A <- insertion_pass(A,M)
        insertion_pass(&A);
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
    std::cout << "Topological pass with _A of size: " <<  _A->size() <<std::endl;
    bool changed = false;
    // E <- edges of the triangles in A
    std::vector<OpenMesh::SmartEdgeHandle> E;
    std::vector<OpenMesh::SmartFaceHandle> newTriangles;

    auto tempMesh = mesh_;
    auto tempA = _A;

    for(unsigned long i=0; i < tempA->size(); ++i){
        auto top = tempA->top();
        tempA->pop();
        auto fh = top.face_handle_;

        for (auto fe_it = tempMesh.fe_cwiter(fh); fe_it.is_valid(); ++fe_it) {
                E.emplace_back(*fe_it);
        }
        // foreach e in E
        // if e still exists attempt to contract e and smooth vertex created
        for(auto e: E){
            if(tempMesh.status(fh).deleted())
                break;

            if (tempMesh.is_flip_ok(e)) {
                tempMesh.flip(e);
                newTriangles = {e.h0().face(), e.h1().face()};
                break;
            }
        }
        E.clear();
    }
    // return surviving set in A and the triangles in M altered by the contractions
    for(auto new_fh: newTriangles){
        if (!mesh_.status(new_fh).deleted()) {
            tempA->push(Triangle(new_fh, QualityEvaluation::evaluate(new_fh,mesh_)));
        }
    }
    if(tempA->top().quality_ > _A->top().quality_){
        *_A = *tempA;
        mesh_ = tempMesh;
        changed = true;
        return changed;
    }

    std::cout << "Topological pass reverted or useless" << std::endl;

    return changed;
}
void MainLoop::edge_contraction_pass(PriorityQueue* _A){
    std::cout << "Edge contraction pass with _A of size: " <<  _A->size() << std::endl;
    // E <- edges of the triangles in A
    std::vector<OpenMesh::SmartEdgeHandle> E;
    std::vector<OpenMesh::SmartFaceHandle> newTriangles;

    for(unsigned long i=0; i < _A->size(); ++i){
        auto top = _A->top();
        _A->pop();
        auto fh = top.face_handle_;

        for (auto fe_it = mesh_.fe_cwiter(fh); fe_it.is_valid(); ++fe_it) {
                E.emplace_back(*fe_it);
        }
        // foreach e in E
        // if e still exists attempt to contract e and smooth vertex created
        for(auto e: E){
            if(mesh_.status(fh).deleted())
                break;

            if (mesh_.is_collapse_ok(e.h0())) {
                auto he = e.h0();
                // get the faces altered by contraction
                for(auto he_fh: he.from().faces()){
                    newTriangles.emplace_back(he_fh);
                }
                mesh_.collapse(he);
                smoother_.smooth(he.to());
                break;
            }
        }
        E.clear();
    }
    // return surviving set in A and the triangles in M altered by the contractions
    for(auto new_fh: newTriangles){
        if (!mesh_.status(new_fh).deleted()) {
            _A->push(Triangle(new_fh, QualityEvaluation::evaluate(new_fh,mesh_)));
        }
    }

}



void MainLoop::insertion_pass(PriorityQueue* _A){
    std::cout << "Insertion pass with _A of size: " <<  _A->size()  << std::endl;

    std::vector<OpenMesh::SmartFaceHandle> facesContainingP;
    // foreach triangle in _A
    for(unsigned long i=0; i < _A->size(); ++i){
        // Reset mesh property before pass
        for(auto face: mesh_.faces()){
            mesh_.property(face_visited_, face) = false;
        }
        auto top = _A->top();
        _A->pop();
        auto fh = top.face_handle_;

        Point p;
        Point vertices[3];
        int j = 0;

        for(auto vh: fh.vertices()){
            vertices[j++] = mesh_.point(vh);
        }
        // find the ideal position of new vertex p
        p = (vertices[0] + vertices[1] + vertices[2])/3;
        // Find all existing triangles whose circumscribing circle contains the new point.
        find_faces_with_p(facesContainingP, fh, p);

        for(auto f: facesContainingP){
            mesh_.delete_face(f);
        }

    }

    // Bowyer-Watson algorithm
    // This can be done to find the triangle which contains the new point first. Then the neighbours of this
//    triangle are searched and then their neighbours, etc., until no more neighbours have the
//    new point in their circumscribing circle.
//    3. Delete these triangles; this creates a convex cavity.
//    4. Join the new point to all the vertices on the boundary of the cavity


}

void MainLoop::find_faces_with_p(std::vector<OpenMesh::SmartFaceHandle>& _list, OpenMesh::SmartFaceHandle _fh, const Point _p){
    if(mesh_.property(face_visited_, _fh) || !contains_p(_fh, _p))
        return;
    _list.emplace_back(_fh);
    mesh_.property(face_visited_, _fh) = true;
    for(auto ffh: _fh.faces()){
        find_faces_with_p(_list, ffh, _p);
    }

}
// from https://stackoverflow.com/questions/39984709/how-can-i-check-wether-a-point-is-inside-the-circumcircle-of-3-points
bool MainLoop::contains_p(OpenMesh::SmartFaceHandle _fh, const Point _p){
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
        std::cout << _fh.idx() << " contains p" << std::endl;
        return true;
    }

    return false;
}

void MainLoop::smoothing_pass(PriorityQueue* _A){
    std::cout << "Smoothing pass with _A of size: " <<  _A->size() << std::endl;
    // V <- the vertices in A
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
        smoother_.smooth(v);
        // A <- A U the triangles adjoining v in M
        for(auto vf_it = mesh_.vf_cwiter(v); vf_it.is_valid(); ++vf_it){
            double quality = QualityEvaluation::evaluate(*vf_it, mesh_);
            _A->push(Triangle(*vf_it,quality));
        }
    }
}
