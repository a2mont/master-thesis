#include "Tests.hh"



bool Tests::t_EdgeRemoval()
{
    std::cout << "\033[1;35m------------------------------"
                 "\nEdge removal test\033[0m" << std::endl;
    TetrahedralMesh mesh;
    TetLoop::PriorityQueue queue;
    EdgeHandle edgeToRemove;
    int cellNb = 0;
    double qualityBefore = 0;
    double qualityAfter = 0;
    std::map<int,int> cv = {{0,0}};

    auto v0 = mesh.add_vertex(ACG::Vec3d( 0,  0,  0));
    auto v1 = mesh.add_vertex(ACG::Vec3d(-1,  5,  1));
    auto v2 = mesh.add_vertex(ACG::Vec3d( 0,  5, -1));
    auto v3 = mesh.add_vertex(ACG::Vec3d( 1,  5,  1));
    auto v4 = mesh.add_vertex(ACG::Vec3d( 0, 10,  0));

    mesh.add_cell(v0,v1,v3,v4, true);
    mesh.add_cell(v0,v2,v1,v4, true);
    mesh.add_cell(v4,v2,v3,v0, true);

    for(auto e: mesh.edges()){
        if(!mesh.is_boundary(e)){
            edgeToRemove = e;
            std::cout << "Removing central edge: "<< e << std::endl;
        }
    }

    computeQuality(queue, mesh);
    qualityBefore = queue.top().quality_;
    TetLoop::reset_queue(queue);

    TetLoop loop(mesh, 0.4, cv);
    loop.edgeRemoval(edgeToRemove);

    computeQuality(queue, mesh);
    qualityAfter = queue.top().quality_;

    for(auto ch: mesh.cells()){
        cellNb++;
    }

    if(qualityAfter <= qualityBefore){
        std::cout << "\033[1;31m X Quality should improve after operation!\033[0m" << std::endl;
        return false;
    }
    if(cellNb != 2){
        std::cout << "\033[1;31m X Only 2 cells should remain!\033[0m" << std::endl;
        return false;
    }

    std::cout << "\033[1;32mPassed\033[0m\n"
                 "\033[1;35m------------------------------\033[0m" << std::endl;

    return true;
}

bool Tests::t_FaceRemoval(){
    std::cout << "\033[1;35m------------------------------"
                 "\nFace removal test\033[0m"<< std::endl;

    std::cout << "\033[1;32mPassed\033[0m\n"
                 "\033[1;35m------------------------------\033[0m" << std::endl;
    return true;
}

bool Tests::t_EdgeContraction(){
    std::cout << "\033[1;35m------------------------------"
                 "\nEdge contraction test\033[0m"<< std::endl;
    TetrahedralMesh mesh;
    TetLoop::PriorityQueue queue;
    EdgeHandle edgeToRemove;
    std::vector<CellHandle> tetsAltered;
    auto v0 = mesh.add_vertex(ACG::Vec3d( 0,  0,  0));
    auto v1 = mesh.add_vertex(ACG::Vec3d(-1,  5,  1));
    auto v2 = mesh.add_vertex(ACG::Vec3d( 0,  5, -1));
    auto v3 = mesh.add_vertex(ACG::Vec3d( 1,  5,  1));
    auto v4 = mesh.add_vertex(ACG::Vec3d( 0, 10,  0));

    mesh.add_cell(v0,v1,v3,v4, true);
    mesh.add_cell(v0,v2,v1,v4, true);
    mesh.add_cell(v4,v2,v3,v0, true);

    edgeToRemove = EdgeHandle(3);
    std::cout << "Edge to contract: "<< edgeToRemove << std::endl;

    std::map<int,int> cv = {{0,0}};

    computeQuality(queue, mesh);

    TetLoop loop(mesh, 0.4, cv);
    loop.contractEdge(edgeToRemove, tetsAltered);
    std::cout << "Tets altered" << std::endl;
    for(auto t: tetsAltered){
        std::cout << t << ", ";
    }

    std::cout << "\033[1;32mPassed\033[0m\n"
                 "\033[1;35m------------------------------\033[0m" << std::endl;

    return true;
}

bool Tests::runAll(){
    std::cout << "\033[1;35mRun all tests\n" << std::endl;
    bool passed, contract, edge, face = false;
    contract    = t_EdgeContraction();
    edge        = t_EdgeRemoval();
    face        = t_FaceRemoval();
    passed = contract && edge && face;
    return passed;
}

void Tests::computeQuality(TetLoop::PriorityQueue& _queue, TetrahedralMesh& _mesh){
    TetLoop::reset_queue(_queue);
    for(auto c_it = _mesh.cells_begin(); c_it != _mesh.cells_end(); ++c_it){
        double quality = QualityEvaluation::evaluate(*c_it, _mesh);
        _queue.push(TetLoop::Tet(*c_it, quality));
    }
}

