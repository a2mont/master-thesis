#include "Smoothing.hh"

void Smoothing::smooth(TriMesh& _mesh, OpenMesh::SmartVertexHandle _vh){
    ACG::Vec3d pt(0.,0.,0.);
    unsigned int n_neighbours = 0;

    if(_vh.is_boundary())
        return;
    for(auto vvh: _vh.vertices()){
        pt += _mesh.point(vvh);
        n_neighbours++;
    }
    pt /= n_neighbours;
    _mesh.point(_vh) = pt;
}

//void Smoothing::smooth(TetrahedralMesh& _mesh, OpenVolumeMesh::VertexHandle _vh){
//    ACG::Vec3d pt(0.,0.,0.);
//    unsigned int n_neighbours = 0;

//    if(_mesh.is_boundary(_vh))
//        return;
//    for(auto vvh: _mesh.vertex_vertices(_vh)){
//        pt += _mesh.vertex(vvh);
//        n_neighbours++;
//    }
//    pt /= n_neighbours;
//    _mesh.set_vertex(_vh, pt);
//}
using namespace OpenVolumeMesh;
void Smoothing::smooth(TetrahedralMesh& _mesh, OpenVolumeMesh::VertexHandle _vh){
    if(_mesh.is_boundary(_vh))
        return;
    bool printDebug(false);
    ACG::Vec3d pt(0.,0.,0.);
    unsigned int n_neighbours = 0;
    std::set<HalfFaceHandle> neighbors_opposite_hf;
    for(auto ch: _mesh.vertex_cells(_vh)){
        HalfFaceHandle opposite_hf;
        for(auto hf: _mesh.cell_halffaces(ch)){
            auto verts = _mesh.get_halfface_vertices(hf);
            if(std::find(verts.begin(), verts.end(), _vh) == verts.end()){
                opposite_hf = hf;
                break;
            }
        }
        auto opp = _mesh.opposite_halfface_handle(opposite_hf);
        if(_mesh.is_valid(opp)){
            neighbors_opposite_hf.insert(opp);
        }
    }

    // 0.1 MUST match CHEBY_TRESHOLD from TetLoop
    bool success = find_chebyshev_center(_mesh, neighbors_opposite_hf, 0.1, pt);
    if(success){
        if(printDebug){
            std::cout << "Smoothing vertex from " << _mesh.vertex(_vh) << " to " << pt << std::endl;
        }
        _mesh.set_vertex(_vh, pt);
    }
}


void Smoothing::triangle_normal_and_centroid(const TetrahedralMesh& _mesh,
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

bool Smoothing::find_chebyshev_center(const TetrahedralMesh& mesh,
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
