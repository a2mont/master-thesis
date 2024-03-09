#include "Tests.hh"

bool Tests::t_EdgeRemoval()
{
    std::cout << "\033[1;35m------------------------------"
                 "\nEdge removal test\n\033[0m" << std::endl;
    TetrahedralMesh mesh;
    TetLoop::PriorityQueue queue;
    EdgeHandle edgeToRemove;
    std::vector<CellHandle> added;
    int cellNb = 0;
    double qualityBefore = 0;
    double qualityAfter = 0;
    double qualityDelta  = 0;
    double rejectedTotal = 0;
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

    TetLoop loop(mesh,0.4, cv, true, false);
    auto begin = std::chrono::high_resolution_clock::now();
    loop.edgeRemoval(edgeToRemove, added, qualityDelta, rejectedTotal, true);
    TetLoop::displayIterationTime(begin, "Edge removal test");


    TetLoop::computeQuality(queue, mesh);
    qualityAfter = queue.top().quality_;

    for(auto ch: mesh.cells()){
        cellNb++;
    }

    std::cout << "Quality: "<< qualityBefore << " -> " << qualityAfter << std::endl;
    std::cout << "\t-Delta: "<< qualityDelta << std::endl;
    if(qualityAfter <= qualityBefore){
        std::cout << "\033[1;31m X Quality should improve after operation!\n\033[0m" << std::endl;
        return false;
    }
    std::cout << "Cells #: "<< cellNb << std::endl;
    if(cellNb != 2){
        std::cout << "\033[1;31m X Only 2 cells should remain!\n\033[0m" << std::endl;
        return false;
    }
    std::cout << "Modified tets: "<< added.size() << std::endl;
    if(added.size() != 2){
        std::cout << "\033[1;31m X The 2 modified tets should be added!\n\033[0m" << std::endl;
        return false;
    }

    std::cout << "\033[1;32mPassed\033[0m\n"
                 "\033[1;35m------------------------------\n\033[0m" << std::endl;

    return true;
}
bool Tests::t_custom_EdgeRemoval(std::string _filename)
{
    std::cout << "\033[1;35m------------------------------"
                 "\nEdge removal test with custom mesh\n\033[0m" << std::endl;
    TetrahedralMesh mesh;
    TetLoop::PriorityQueue queue;
    EdgeHandle edgeToRemove;
    std::vector<CellHandle> added;
    int cellNb = 0;
    double qualityBefore = 0;
    double qualityAfter  = 0;
    double qualityDelta  = 0;
    double rejectedTotal = 0;
    std::map<int,int> cv = {{0,0}};

    OpenVolumeMesh::IO::FileManager fileManager;
    fileManager.readFile(LOGS_MESH + _filename, mesh);

    for(auto e: mesh.edges()){
        if(!mesh.is_boundary(e)){
            edgeToRemove = e;
            std::cout << "Removing central edge: "<< e << std::endl;
        }
    }

    TetLoop::computeQuality(queue, mesh);
    qualityBefore = queue.top().quality_;

    TetLoop loop(mesh, 0.4, cv, true, false);
    auto begin = std::chrono::high_resolution_clock::now();
    loop.edgeRemoval(edgeToRemove, added, qualityDelta, rejectedTotal, true);
    TetLoop::displayIterationTime(begin, "Edge removal test");

    TetLoop::computeQuality(queue, mesh);
    qualityAfter = queue.top().quality_;

    for(auto ch: mesh.cells()){
        cellNb++;
    }

    std::cout << "Quality: "<< qualityBefore << " -> " << qualityAfter << std::endl;
    std::cout << "\t-Delta: "<< qualityDelta << std::endl;
    if(qualityAfter <= qualityBefore){
        std::cout << "\033[1;31m X Quality should improve after operation!\n\033[0m" << std::endl;
        return false;
    }
    std::cout << "Cells #: "<< cellNb << std::endl;
    if(cellNb != 2){
        std::cout << "\033[1;31m X Only 2 cells should remain!\n\033[0m" << std::endl;
        return false;
    }
    std::cout << "Modified tets: "<< added.size() << std::endl;
    if(added.empty()){
        std::cout << "\033[1;31m X The modified tets should be added!\n\033[0m" << std::endl;
        return false;
    }

    std::cout << "\033[1;32mPassed\033[0m\n"
                 "\033[1;35m------------------------------\n\033[0m" << std::endl;

    return true;
}

bool Tests::t_FaceRemoval(){
    std::cout << "\033[1;35m------------------------------"
                 "\nFace removal test\n\033[0m"<< std::endl;
    TetrahedralMesh mesh;
    TetLoop::PriorityQueue queue;
    int cellNb = 0;
    double qualityBefore = 0;
    double qualityAfter = 0;
    double qualityDelta  = 0;
    double rejectedTotal = 0;
    std::map<int,int> cv = {{0,0}};

    auto v0 = mesh.add_vertex(ACG::Vec3d( 0,  -10,  0));
    auto v1 = mesh.add_vertex(ACG::Vec3d(-10,  0.5,  10));
    auto v2 = mesh.add_vertex(ACG::Vec3d( 0,  0.5, -10));
    auto v3 = mesh.add_vertex(ACG::Vec3d( 10,  0.5,  10));
    auto v4 = mesh.add_vertex(ACG::Vec3d( 0, 1,  0));


    mesh.add_cell(v1,v2,v3,v0, true);
    mesh.add_cell(v1,v3,v2,v4, true);

    TetLoop::computeQuality(queue, mesh);
    qualityBefore = queue.top().quality_;
    TetLoop::reset_queue(queue);

    TetLoop loop(mesh, 0.4, cv, true, false);
    auto begin = std::chrono::high_resolution_clock::now();
    loop.faceRemoval(FaceHandle(0), qualityDelta, rejectedTotal, false);
    TetLoop::displayIterationTime(begin, "Face removal test");

    TetLoop::computeQuality(queue, mesh);
    qualityAfter = queue.top().quality_;

    for(auto ch: mesh.cells()){
        cellNb++;
    }

    std::cout << "Quality: "<< qualityBefore << " -> " << qualityAfter << std::endl;
    if(qualityAfter <= qualityBefore){
        std::cout << "\033[1;31m X Quality should improve after operation!\n\033[0m" << std::endl;
        return false;
    }
    if(cellNb != 3){
        std::cout << "\033[1;31m X Only 3 cells should remain!\n\033[0m" << std::endl;
        return false;
    }

    std::cout << "\033[1;32mPassed\033[0m\n"
                 "\033[1;35m------------------------------\n\033[0m" << std::endl;
    return true;
}

bool Tests::t_flip23(){
    std::cout << "\033[1;35m------------------------------"
                 "\n2-3flip test\n\033[0m"<< std::endl;
    TetrahedralMesh mesh;
    TetLoop::PriorityQueue queue;
    FaceHandle toRemove(-1);
    int cellNb = 0;
    double qualityBefore = 0;
    double qualityAfter = 0;
    double qualityDelta  = 0;
    std::map<int,int> cv = {{0,0}};

    auto v0 = mesh.add_vertex(ACG::Vec3d( 0,  -10,  0));
    auto v1 = mesh.add_vertex(ACG::Vec3d(-10,  0.5,  10));
    auto v2 = mesh.add_vertex(ACG::Vec3d( 0,  0.5, -10));
    auto v3 = mesh.add_vertex(ACG::Vec3d( 10,  0.5,  10));
    auto v4 = mesh.add_vertex(ACG::Vec3d( 0, 1,  0));


    mesh.add_cell(v1,v2,v3,v0, true);
    mesh.add_cell(v1,v3,v2,v4, true);

    TetLoop::computeQuality(queue, mesh);
    qualityBefore = queue.top().quality_;
    TetLoop::reset_queue(queue);

    for(auto fh: mesh.faces()){
        if(!mesh.is_boundary(fh)){
            std::cout << "Face to remove: "<< fh << std::endl;
            toRemove = fh;
            break;
        }
    }

    std::vector<CellHandle> added;
    TetLoop loop(mesh, 0.4, cv, true, false);
    auto begin = std::chrono::high_resolution_clock::now();
    loop.flip23(mesh, toRemove, added);
    TetLoop::displayIterationTime(begin, "2-3 flip test");

    TetLoop::computeQuality(queue, mesh);
    qualityAfter = queue.top().quality_;

    for(auto ch: mesh.cells()){
        cellNb++;
    }

    std::cout << "Quality: "<< qualityBefore << " -> " << qualityAfter << std::endl;
    if(qualityAfter <= qualityBefore){
        std::cout << "\033[1;31m X Quality should improve after operation! Quality: "
                  << qualityBefore<< " -> "<< qualityAfter <<"\n\033[0m" << std::endl;
        return false;
    }
    if(cellNb != 3){
        std::cout << "\033[1;31m X Only 3 cells should remain! Cell number:"<< cellNb
                  << "\n\033[0m" << std::endl;
        return false;
    }

    std::cout << "Added cells: "<< TetLoop::vectorToString<CellHandle>(added) << std::endl;
    if(added.size() != 3){
        std::cout << "\033[1;31m X Only 3 cells should be added ! added: "<< added.size()
                  <<"\n\033[0m" << std::endl;
        return false;
    }

    std::cout << "\033[1;32mPassed\033[0m\n"
                 "\033[1;35m------------------------------\n\033[0m" << std::endl;
    return true;
}

bool Tests::t_flip32(){
    std::cout << "\033[1;35m------------------------------"
                 "\n3-2 flip test\n\033[0m" << std::endl;
    TetrahedralMesh mesh;
    TetLoop::PriorityQueue queue;
    EdgeHandle edgeToRemove;
    std::vector<CellHandle> added;
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

    TetLoop loop(mesh, 0.4, cv, true, false);
    auto begin = std::chrono::high_resolution_clock::now();
    loop.flip32(mesh,edgeToRemove, added);
    TetLoop::displayIterationTime(begin, "3-2 flip test");


    TetLoop::computeQuality(queue, mesh);
    qualityAfter = queue.top().quality_;

    for(auto ch: mesh.cells()){
        cellNb++;
    }

    std::cout << "Quality: "<< qualityBefore << " -> " << qualityAfter << std::endl;
    if(qualityAfter <= qualityBefore){
        std::cout << "\033[1;31m X Quality should improve after operation! Quality: "
                  << qualityBefore<< " -> "<< qualityAfter <<"\n\033[0m" << std::endl;
        return false;
    }
    if(cellNb != 2){
        std::cout << "\033[1;31m X Only 2 cells should remain! Cell number:"<< cellNb
                  << "\n\033[0m" << std::endl;
        return false;
    }

    std::cout << "Added cells: "<< TetLoop::vectorToString<CellHandle>(added) << std::endl;
    if(added.size() != 2){
        std::cout << "\033[1;31m X Only 2 cells should be added ! added: "<< added.size()
                  <<"\n\033[0m" << std::endl;
        return false;
    }

    std::cout << "\033[1;32mPassed\033[0m\n"
                 "\033[1;35m------------------------------\n\033[0m" << std::endl;

    return true;
}

bool Tests::t_multiface(){
    std::cout << "\033[1;35m------------------------------"
                 "\nMultiface removal test\n\033[0m"<< std::endl;
    TetrahedralMesh mesh;
    TetLoop::PriorityQueue queue;
    FaceHandle toRemove(-1);
    int cellNb = 0;
    double qualityBefore = 0;
    double qualityAfter = 0;
    std::map<int,int> cv = {{0,0}};

    auto v0 = mesh.add_vertex(ACG::Vec3d( 0,  -10,  0));
    auto v1 = mesh.add_vertex(ACG::Vec3d(-10,  0.5,  10));
    auto v2 = mesh.add_vertex(ACG::Vec3d( 0,  0.5, -10));
    auto v3 = mesh.add_vertex(ACG::Vec3d( 10,  0.5,  10));
    auto v4 = mesh.add_vertex(ACG::Vec3d( 0, 1,  0));


    mesh.add_cell(v1,v2,v3,v0, true);
    mesh.add_cell(v1,v3,v2,v4, true);

//    double quality = QualityEvaluation::evaluate(v0, v4, v3, v2, mesh);
//    std::cout << "Quality: "<< quality << std::endl;
//    quality = QualityEvaluation::evaluate(v0, v4, v1, v3, mesh);
//    std::cout << "Quality: "<< quality << std::endl;
//    quality = QualityEvaluation::evaluate(v0, v4, v2, v1, mesh);
//    std::cout << "Quality: "<< quality << std::endl;

    TetLoop::computeQuality(queue, mesh);
    qualityBefore = queue.top().quality_;
    TetLoop::reset_queue(queue);

    for(auto fh: mesh.faces()){
        if(!mesh.is_boundary(fh)){
            std::cout << "Face to remove: "<< fh << std::endl;
            toRemove = fh;
            break;
        }
    }

    std::vector<CellHandle> added;
    TetLoop loop(mesh, 0.4, cv, true, false);
    auto begin = std::chrono::high_resolution_clock::now();
    loop.multiFace(mesh, toRemove, added);
    TetLoop::displayIterationTime(begin, "Multiface flip test");


    TetLoop::computeQuality(queue, mesh);
    qualityAfter = queue.top().quality_;

    for(auto ch: mesh.cells()){
        cellNb++;
    }

    std::cout << "Quality: "<< qualityBefore << " -> " << qualityAfter << std::endl;
    if(qualityAfter <= qualityBefore){
        std::cout << "\033[1;31m X Quality should improve after operation! Quality: "
                  << qualityBefore<< " -> "<< qualityAfter <<"\n\033[0m" << std::endl;
        return false;
    }
    if(cellNb != 3){
        std::cout << "\033[1;31m X Only 3 cells should remain! Cell number:"<< cellNb
                  << "\n\033[0m" << std::endl;
        return false;
    }

    std::cout << "Added cells: "<< TetLoop::vectorToString<CellHandle>(added) << std::endl;
    if(added.size() != 3){
        std::cout << "\033[1;31m X Only 3 cells should be added ! added: "<< added.size()
                  <<"\n\033[0m" << std::endl;
        return false;
    }

    std::cout << "\033[1;32mPassed\033[0m\n"
                 "\033[1;35m------------------------------\n\033[0m" << std::endl;
    return true;
}

bool Tests::t_custom_multiface(std::string _filename){
    std::cout << "\033[1;35m------------------------------"
                 "\nMultiface removal test\n\033[0m"<< std::endl;
    TetrahedralMesh mesh;

    OpenVolumeMesh::IO::FileManager fileManager;
    fileManager.readFile(LOGS_MESH + _filename, mesh);

    TetLoop::PriorityQueue queue;
    FaceHandle toRemove(-1);
    int cellNb = 0;
    double qualityBefore = 0;
    double qualityAfter = 0;
    std::map<int,int> cv = {{0,0}};

//    double quality = QualityEvaluation::evaluate(v0, v4, v3, v2, mesh);
//    std::cout << "Quality: "<< quality << std::endl;
//    quality = QualityEvaluation::evaluate(v0, v4, v1, v3, mesh);
//    std::cout << "Quality: "<< quality << std::endl;
//    quality = QualityEvaluation::evaluate(v0, v4, v2, v1, mesh);
//    std::cout << "Quality: "<< quality << std::endl;

    TetLoop::computeQuality(queue, mesh);
    qualityBefore = queue.top().quality_;
    TetLoop::reset_queue(queue);

    for(auto fh: mesh.faces()){
        if(!mesh.is_boundary(fh)){
            std::cout << "Face to remove: "<< fh << std::endl;
            toRemove = fh;
            break;
        }
    }

    std::vector<CellHandle> added;
    TetLoop loop(mesh, 0.4, cv);
    auto begin = std::chrono::high_resolution_clock::now();
    loop.multiFace(mesh, toRemove, added);
    TetLoop::displayIterationTime(begin, "Multiface flip test");


    TetLoop::computeQuality(queue, mesh);
    qualityAfter = queue.top().quality_;

    for(auto ch: mesh.cells()){
        cellNb++;
    }

    std::cout << "Quality: "<< qualityBefore << " -> " << qualityAfter << std::endl;
    if(qualityAfter <= qualityBefore){
        std::cout << "\033[1;31m X Quality should improve after operation! Quality: "
                  << qualityBefore<< " -> "<< qualityAfter <<"\n\033[0m" << std::endl;
        return false;
    }
    if(cellNb != 3){
        std::cout << "\033[1;31m X Only 3 cells should remain! Cell number:"<< cellNb
                  << "\n\033[0m" << std::endl;
        return false;
    }

    std::cout << "Added cells: "<< TetLoop::vectorToString<CellHandle>(added) << std::endl;
    if(added.size() != 3){
        std::cout << "\033[1;31m X Only 3 cells should be added ! added: "<< added.size()
                  <<"\n\033[0m" << std::endl;
        return false;
    }

    std::cout << "\033[1;32mPassed\033[0m\n"
                 "\033[1;35m------------------------------\n\033[0m" << std::endl;
    return true;
}

bool Tests::t_EdgeContraction(){
    std::cout << "\033[1;35m------------------------------"
                 "\nEdge contraction test\n\033[0m"<< std::endl;
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
    double delta = 0;
    std::cout << "Tets before: " << mesh.n_logical_cells() << std::endl;
    loop.contractEdge(edgeToRemove, tetsAltered, delta, delta);
    std::cout << "Tets altered" << std::endl;
    for(auto t: tetsAltered){
        std::cout << t << ", ";
    }
    std::cout << "Remaining tets: " << mesh.n_logical_cells() << std::endl;

    std::cout << "\033[1;32mPassed\033[0m\n"
                 "\033[1;35m------------------------------\n\033[0m" << std::endl;

    return true;
}

bool Tests::t_chebyshev_centroid(){
    std::cout << "\033[1;35m------------------------------"
                 "\nChebyshev center test\n\033[0m"<< std::endl;
    TetrahedralMesh mesh;

    auto v0 = mesh.add_vertex({1.5, 1.5, -1}); //cluster virtual vertex
    auto v1 = mesh.add_vertex({0,2,0});
    auto v2 = mesh.add_vertex({1,0,0});
    auto v3 = mesh.add_vertex({1,1,0});
    auto v4 = mesh.add_vertex({2,1,1});
    auto v5 = mesh.add_vertex({2,0,0});
    auto v6 = mesh.add_vertex({3,2,0});


    mesh.add_cell({v3, v2, v1, v0});
    mesh.add_cell({v3, v1, v6, v0});
    mesh.add_cell({v3, v6, v4, v0});
    mesh.add_cell({v4, v6, v5, v0});


    std::map<int,int> cv = {{0,0}};
    TetLoop loop(mesh, 0.4, cv);
    std::set<HalfFaceHandle> hfs;
    for(auto hf:mesh.halffaces()){
        if(mesh.is_boundary(hf)){
            hfs.insert(hf);
        }
    }

    ACG::Vec3d new_pos;
    auto begin = std::chrono::high_resolution_clock::now();
    bool cheb_result =
            loop.find_chebyshev_center(mesh, hfs, TetLoop::CHEBY_THRESHOLD, new_pos);
    TetLoop::displayIterationTime(begin, "Cheby center test");


    if(!cheb_result){
        std::cout << "\033[1;31m X Failed to find cheby center \n\033[0m" << std::endl;
        return false;
    }

    std::cout<<" --> new position = "<<new_pos<<std::endl;
    std::cout << "\033[1;32mPassed\033[0m\n"
                 "\033[1;35m------------------------------\n\033[0m" << std::endl;
    return true;
}

bool Tests::t_custom_chebyshev_centroid(std::string _filename){
    std::cout << "\033[1;35m------------------------------"
                 "\nChebyshev center test on file: "<< _filename <<"\n\033[0m"<< std::endl;
    TetrahedralMesh mesh;

    OpenVolumeMesh::IO::FileManager fileManager;
    fileManager.readFile(LOGS_MESH + _filename, mesh);

    std::map<int,int> cv = {{0,0}};
    TetLoop loop(mesh, 0.5, cv);
    std::set<HalfFaceHandle> hfs;
    for(auto hf:mesh.halffaces()){
        if(mesh.is_boundary(hf)){
            hfs.insert(hf);
        }
    }
    std::cout << "Bounds size:" << hfs.size()<< std::endl;


    ACG::Vec3d new_pos;
    auto cheb_result = loop.find_chebyshev_center(mesh, hfs, TetLoop::CHEBY_THRESHOLD, new_pos);

    if(!cheb_result){
        std::cout<<" --> failed to find Chebyshev center"<<std::endl;
        return false;
    }

    std::cout<<" --> new position = "<<new_pos<<std::endl;
    std::cout << "\033[1;32mPassed\033[0m\n"
                 "\033[1;35m------------------------------\n\033[0m" << std::endl;
    return true;
}


bool Tests::t_quality_evaluation(){
    std::cout << "\033[1;35m------------------------------"
                 "\nQuality evalutation test\n\033[0m"<< std::endl;
    TetrahedralMesh mesh;
    TetLoop::PriorityQueue queue;
    // Good setting
    auto v0 = mesh.add_vertex({0,0,0});
    auto v1 = mesh.add_vertex({0,1,1});
    auto v2 = mesh.add_vertex({1,0,1});
    auto v3 = mesh.add_vertex({1,1,0});
    mesh.add_cell({v0,v1,v2,v3});

    TetLoop::computeQuality(queue,mesh);
    double qualityBefore = queue.top().quality_;
    std::cout << "Quality of good element: " << qualityBefore <<std::endl;

    mesh.clear();

    // Bad setting
    v0 = mesh.add_vertex({0,0,0});
    v1 = mesh.add_vertex({0,1,10});
    v2 = mesh.add_vertex({1,0,0.1});
    v3 = mesh.add_vertex({1,5,0});
    mesh.add_cell({v0,v1,v2,v3});


    TetLoop::computeQuality(queue,mesh);
    double qualityAfter = queue.top().quality_;
    std::cout << "Quality of bad element: "<< qualityAfter << std::endl;

    if(qualityAfter > qualityBefore){
        std::cout << "\033[1;31m X Quality should improve!\n\033[0m" << std::endl;
        return false;
    }

    mesh.clear();
    TetrahedralizedVoxelGridGenerator<TetrahedralMesh>::
            generate_mesh(3,3,3,mesh);
    std::map<int,int> cst;
    Experiment3D exp(mesh, 0.1, cst, false);
    exp.generate_torsion_mesh(1.5/15 *5, false);

    TetLoop::computeQuality(queue,mesh);
    int counter(0);
    std::cout << "Ordered list of elements: " << std::endl;
    while (!queue.empty()){
        auto top = queue.top();
        queue.pop();
        std::cout << top.quality_ << ", ";
        if(++counter % 10 == 0){
            std::cout << std::endl;
        }
    }
    std::cout << std::endl;

    std::cout << "\033[1;32mPassed\033[0m\n"
                 "\033[1;35m------------------------------\n\033[0m" << std::endl;

    return true;
}

bool Tests::t_speed(){
    std::cout << "\033[1;35m------------------------------"
                 "\nSpeed test\n\033[0m"<< std::endl;
    auto begin = std::chrono::high_resolution_clock::now();
    TetrahedralMesh mesh;
    auto mesh_begin = std::chrono::high_resolution_clock::now();
    TetrahedralizedVoxelGridGenerator<TetrahedralMesh>::
            generate_mesh(10, 10,15 , mesh);
    TetLoop::displayIterationTime(mesh_begin, "10x10x15 Mesh generation");

    TetrahedralMesh temp;
    auto copy_begin = std::chrono::high_resolution_clock::now();
    temp = mesh;
    TetLoop::displayIterationTime(copy_begin, "Mesh copy");

    TetLoop::PriorityQueue queue;
    auto quality_begin = std::chrono::high_resolution_clock::now();
    TetLoop::computeQuality(queue, mesh);
    TetLoop::displayIterationTime(quality_begin, "Whole mesh quality evaluation");

    TetLoop::reset_queue(queue);
    std::set<CellHandle> cellSet({CellHandle(0), CellHandle(10), CellHandle(105),
                                  CellHandle(666), CellHandle(354), CellHandle(21),
                                  CellHandle(78), CellHandle(922), CellHandle(545),
                                  CellHandle(123), CellHandle(330), CellHandle(677)});
    auto quality_set_begin = std::chrono::high_resolution_clock::now();
    for(auto ch: cellSet){
        // minus to handle the change of quality to dirichlet
        double quality = QualityEvaluation::evaluate(ch, mesh);
        queue.push(TetLoop::Tet(ch, quality));
    }
    TetLoop::displayIterationTime(quality_set_begin, "Small Set quality evaluation");

    TetLoop::reset_queue(queue);
    std::set<CellHandle> bigCellSet;
    for(int i = 0; i < 1500; ++i){
        bigCellSet.insert(CellHandle(i));
    }
    quality_set_begin = std::chrono::high_resolution_clock::now();
    for(auto ch: bigCellSet){
        // minus to handle the change of quality to dirichlet
        double quality = QualityEvaluation::evaluate(ch, mesh);
        queue.push(TetLoop::Tet(ch, quality));
    }
    TetLoop::displayIterationTime(quality_set_begin, "Big Set quality evaluation");

    TetLoop::displayIterationTime(begin, "Complete test");
    std::cout << "\033[1;32mDone\033[0m\n"
                 "\033[1;35m------------------------------\n\033[0m" << std::endl;

    return true;
}

bool Tests::runAll(){
    std::cout << "\033[1;35mRun all tests\n" << std::endl;
    bool passed(false), edge(false), multiface(false),
            flip23(false), flip32(false),cheby(false),speed(false);

    edge        = t_EdgeRemoval();
    flip23      = t_flip23();
    flip32      = t_flip32();
    multiface   = t_multiface();
    cheby       = t_chebyshev_centroid();
    speed       = t_speed();
    passed =
            edge        &&
            multiface   &&
            cheby       &&
            flip32      &&
            flip23      &&
            speed
            ;
    return passed;
}


void Tests::errorReproduction(){
    TetrahedralMesh baseMesh;
    auto v0 = baseMesh.add_vertex(ACG::Vec3d( 0,  0,  0));
    auto v1 = baseMesh.add_vertex(ACG::Vec3d(-1,  5,  1));
    auto v2 = baseMesh.add_vertex(ACG::Vec3d( 0,  5, -1));
    auto v3 = baseMesh.add_vertex(ACG::Vec3d( 1,  5,  1));
    auto v4 = baseMesh.add_vertex(ACG::Vec3d( 0, 10,  0));

    baseMesh.add_cell(v0,v1,v3,v4, true);
    baseMesh.add_cell(v0,v2,v1,v4, true);
    baseMesh.add_cell(v4,v2,v3,v0, true);

    TetrahedralMesh& mesh = baseMesh;
    TetrahedralMesh meshCopy = baseMesh;
    TetrahedralMesh copyOfCopy = meshCopy;

    EdgeHandle toSplit(4);
    HalfEdgeHandle toCollapse(22);

    std::cout << "---- Should not work ----" << std::endl;
    copyOfCopy.split_edge(toSplit);
    copyOfCopy.collapse_edge(toCollapse);
    for(auto ch: copyOfCopy.cells()){
        copyOfCopy.get_cell_vertices(ch);
    }

    std::cout << "---- Should work -----" << std::endl;
    mesh.split_edge(toSplit);
    mesh.collapse_edge(toCollapse);
    for(auto ch: mesh.cells()){
        mesh.get_cell_vertices(ch);
    }







}

