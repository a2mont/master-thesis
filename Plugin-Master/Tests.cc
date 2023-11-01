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

    TetLoop::computeQuality(queue, mesh);
    qualityBefore = queue.top().quality_;
    TetLoop::reset_queue(queue);

    TetLoop loop(mesh, 0.4, cv);
    loop.edgeRemoval(edgeToRemove);

    TetLoop::computeQuality(queue, mesh);
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

bool Tests::t_StressEdgeRemoval(){
    std::cout << "\033[1;35m------------------------------"
                 "\nEdge removal stress test\033[0m" << std::endl;

    TetrahedralMesh mesh;
    std::map<int,int> cv = {{0,0}};
    int reverts = 0;
    int iters = 0;
    TetrahedralizedVoxelGridGenerator<TetrahedralMesh>::generate_mesh(5, mesh);

    TetLoop loop(mesh, 0.4, cv);
    for(auto e: mesh.edges()){
        if(!mesh.is_boundary(e)){
            auto result = loop.edgeRemoval(e, false);
            reverts = result ? reverts : reverts + 1;
        }
        if(iters++ % 20 == 0){
            std::cout << "Iter " << iters << "/" << mesh.n_logical_edges() << std::endl;
        }
    }
    std::cout << "# Reverts: " << reverts << "/" << iters << std::endl;
    std::cout << "\033[1;32mPassed\033[0m\n"
                 "\033[1;35m------------------------------\033[0m" << std::endl;

    return true;
}

bool Tests::t_FaceRemoval(){
    std::cout << "\033[1;35m------------------------------"
                 "\nFace removal test\033[0m"<< std::endl;
    TetrahedralMesh mesh;
    TetLoop::PriorityQueue queue;
    int cellNb = 0;
    double qualityBefore = 0;
    double qualityAfter = 0;
    std::map<int,int> cv = {{0,0}};

    auto v0 = mesh.add_vertex(ACG::Vec3d( 0,  0,  0));
    auto v1 = mesh.add_vertex(ACG::Vec3d(-1,  5,  1));
    auto v2 = mesh.add_vertex(ACG::Vec3d( 0,  5, -1));
    auto v3 = mesh.add_vertex(ACG::Vec3d( 1,  5,  1));
    auto v4 = mesh.add_vertex(ACG::Vec3d( 0, 10,  0));


    mesh.add_cell(v1,v2,v3,v0, true);
    mesh.add_cell(v1,v3,v2,v4, true);

    TetLoop::computeQuality(queue, mesh);
    qualityBefore = queue.top().quality_;
    TetLoop::reset_queue(queue);

    TetLoop loop(mesh, 0.4, cv);
    loop.faceRemoval(FaceHandle(0));


    TetLoop::computeQuality(queue, mesh);
    qualityAfter = queue.top().quality_;

    for(auto ch: mesh.cells()){
        cellNb++;
    }

    if(qualityAfter <= qualityBefore){
        std::cout << "\033[1;31m X Quality should improve after operation!\033[0m" << std::endl;
        return false;
    }
    if(cellNb != 3){
        std::cout << "\033[1;31m X Only 2 cells should remain!\033[0m" << std::endl;
        return false;
    }

    std::cout << "\033[1;32mPassed\033[0m\n"
                 "\033[1;35m------------------------------\033[0m" << std::endl;
    return true;
}

bool Tests::t_StressFaceRemoval(){
    std::cout << "\033[1;35m------------------------------"
                 "\nEdge removal stress test\033[0m" << std::endl;

    TetrahedralMesh mesh;
    std::map<int,int> cv = {{0,0}};
    int reverts = 0;
    int iters = 0;
    TetLoop::PriorityQueue queue;
    TetrahedralizedVoxelGridGenerator<TetrahedralMesh>::generate_mesh(5, mesh);

    TetLoop loop(mesh, 0.4, cv);
    int max = mesh.n_logical_faces();
    for(int i = 0; i < max; ++i){
        auto f = FaceHandle(i);
        if(!mesh.is_deleted(f)){
            auto result = loop.faceRemoval(f, false);
            reverts = result ? reverts : reverts + 1;
        }
        if(iters++ % 50 == 0){
            std::cout << "Iter " << iters << "/" << max << std::endl;
        }
    }

    TetLoop::computeQuality(queue, mesh);
    if(queue.top().quality_ <= 10e-5){
        std::cout << "\033[1;31m X creates inverted tets!\033[0m" << std::endl;
        return false;
    }
    std::cout << "# Reverts: " << reverts << "/" << iters << std::endl;
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

    TetLoop::computeQuality(queue, mesh);

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
    bool passed, contract, edge, face, stressEdge, stressFace;
    passed = contract = edge = face = stressEdge = stressFace = false;

    contract    = t_EdgeContraction();
    edge        = t_EdgeRemoval();
//    face        = t_FaceRemoval();
    stressEdge  = t_StressEdgeRemoval();
    stressFace  = t_StressFaceRemoval();
    passed =
            contract &&
            edge &&
            face &&
            stressEdge &&
            stressFace
            ;
    return passed;
}


