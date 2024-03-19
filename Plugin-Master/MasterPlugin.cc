#include "MasterPlugin.hh"

#include <ObjectTypes/PolyMesh/PolyMesh.hh>
#include <ObjectTypes/TriangleMesh/TriangleMesh.hh>
#include "OpenFlipper/BasePlugin/PluginFunctions.hh"


void MasterPlugin::initializePlugin()
{  
    tool_ = new MasterToolbar();
    QSize size(600, 300);
    tool_->resize(size);

    // Connections
    connect(tool_->generateButton, SIGNAL(clicked()), this, SLOT(slot_generate_base_mesh()));
    connect(tool_->selectButton, SIGNAL(clicked()), this, SLOT(slot_show_constraint_vertex()));
    connect(tool_->displaceButton,  SIGNAL(clicked()), this, SLOT(slot_displace_constraint_vertex()));
    connect(tool_->showQualityCheckbox, SIGNAL(toggled(bool)), this, SLOT(slot_show_quality()));
    connect(tool_->clearButton, SIGNAL(clicked()), this, SLOT(slot_clear_constraints()));
    connect(tool_->beginExpButton, SIGNAL(clicked()), this, SLOT(slot_start_experiment()));
    connect(tool_->nextButton, SIGNAL(clicked()), this, SLOT(slot_experiment_loop()));
    connect(tool_->completeButton, SIGNAL(clicked()), this, SLOT(slot_complete_experiment()));


    // Add the Toolbox
    emit addToolbox("Master", tool_, new QIcon());
}

void MasterPlugin::pluginsInitialized(){
    // Arbitrary id for constraint vertex
    constraint_vhs_[0] = 0;
    slot_generate_base_mesh();
//    Tests::runAll();
//    std::cout << "Volume------------------" << std::endl;
//    for(auto ch: tetmesh_.cells()){
//        double vol = QualityEvaluation::calculate_volume(tetmesh_,ch);
//        if(vol < 1e-5){
//            std::cout << "text: "<< vol << std::endl;
//        }
//    }

    emit addPickMode("Pick Constraint");

}

void MasterPlugin::slotMouseEvent(QMouseEvent* _event) {
    // control modifier is reserved for selecting target
    if (_event->modifiers() & (Qt::ControlModifier)) {
        return;
    }
    if(_event->type() != QEvent::MouseButtonPress) {
        return;
    }

    if (PluginFunctions::actionMode() == Viewer::PickingMode) {
        // handle mouse events
        if (_event->button() == Qt::LeftButton) {
            size_t node_idx, target_idx;
            ACG::Vec3d hit_point;
            // pick vertices
            if (PluginFunctions::scenegraphPick(ACG::SceneGraph::PICK_VERTEX, _event->pos(),
                                                node_idx, target_idx, &hit_point)) {
                BaseObjectData *obj;
                if (PluginFunctions::getPickedObject(node_idx, obj)) {
                    // is picked object a triangle mesh?
                    TriMeshObject *tri_obj = PluginFunctions::triMeshObject(obj);
                    if (tri_obj) {
                        auto vh = tri_obj->mesh()->vertex_handle(target_idx);
                        if (vh == TriMesh::InvalidVertexHandle)
                            return;

                        if(constraint_vhs_.count(vh.idx()) == 0)
                            constraint_vhs_[vh.idx()] = vh.idx();
                        else
                            constraint_vhs_.erase(vh.idx());

                        std::cout << "Constraint vertices" << std::endl;
                        for(auto v: constraint_vhs_){
                            std::cout << v.first <<"," ;
                        }
                        std::cout << "\nTotal: " << constraint_vhs_.size() << std::endl;
                        slot_show_constraint_vertex();

                        return;
                    }
                    // is picked object a tet mesh?
                    TetrahedralMeshObject *tet_obj = PluginFunctions::tetrahedralMeshObject(obj);
                    if (tet_obj) {
                        auto vh = OpenVolumeMesh::VertexHandle(target_idx);
                        if (vh == TetrahedralMesh::InvalidVertexHandle)
                            return;

                        if(constraint_vhs_.count(vh.idx()) == 0)
                            constraint_vhs_[vh.idx()] = vh.idx();
                        else
                            constraint_vhs_.erase(vh.idx());

                        std::cout << "Constraint vertices" << std::endl;
                        for(auto v: constraint_vhs_){
                            std::cout << v.first <<"," ;
                        }
                        std::cout << "\nTotal: " << constraint_vhs_.size() << std::endl;
//                        slot_show_constraint_vertex();

                        return;
                    }
                }
            }
        }
    }

    emit updateView();
}

void MasterPlugin::slot_clear_constraints(){
    constraint_vhs_.clear();
    std::cout << "Constraint vertices cleared" << std::endl;
    slot_show_constraint_vertex();
}

void MasterPlugin::slot_complete_experiment(){
    int remain = timesteps_ - t_;
    std::cout << "Remaining timesteps: "<< remain << std::endl;
    for(remain; remain > 0; --remain){
        tool_->nextButton->click();
    }

}
void MasterPlugin::slot_show_quality(){
    for (PluginFunctions::ObjectIterator o_it(PluginFunctions::TARGET_OBJECTS, DATA_TRIANGLE_MESH);
         o_it != PluginFunctions::objectsEnd(); ++o_it) {
        auto tri_obj = PluginFunctions::triMeshObject(*o_it);
        auto trimesh = tri_obj->mesh();


        if(tool_->showQualityCheckbox->isChecked()){
            Highlight::highlight_triangles(*trimesh);
        }else{
            for(auto fh: trimesh->faces()){
                trimesh->set_color(fh,ACG::Vec4f(1,1,1,1));
            }
        }

        tri_obj->meshNode()->drawMode(
                    ACG::SceneGraph::DrawModes::WIREFRAME
                  | ACG::SceneGraph::DrawModes::POINTS_COLORED
                  | ACG::SceneGraph::DrawModes::SOLID_FACES_COLORED);
        tri_obj->materialNode()->enable_alpha_test(0.8);


        emit updatedObject(tri_obj->id(), UPDATE_COLOR);
    }

}

void MasterPlugin::slot_show_constraint_vertex(){
    if(tool_->selectButton->isChecked()) {
        // Picking mode of our plugin shall be activated
        // set OpenFlipper's action mode to picking
        PluginFunctions::actionMode( Viewer::PickingMode );
        // Now activate our picking mode
        PluginFunctions::pickMode( "Pick constraint" );
    } else {
        // Picking mode shall be deactivated
        PluginFunctions::actionMode( Viewer::ExamineMode );
    }

    for (PluginFunctions::ObjectIterator o_it(PluginFunctions::TARGET_OBJECTS, DATA_TRIANGLE_MESH);
         o_it != PluginFunctions::objectsEnd(); ++o_it) {
        auto tri_obj = PluginFunctions::triMeshObject(*o_it);
        auto trimesh = tri_obj->mesh();

        tri_obj->materialNode()->set_point_size(12);

        highlight_constraints_vertices(trimesh);

        tri_obj->meshNode()->drawMode(
                    ACG::SceneGraph::DrawModes::WIREFRAME
                  | ACG::SceneGraph::DrawModes::POINTS_COLORED
                  | ACG::SceneGraph::DrawModes::SOLID_FACES_COLORED);

        tri_obj->materialNode()->enable_alpha_test(0.8);


        emit updatedObject(tri_obj->id(), UPDATE_COLOR);
    }
    for (PluginFunctions::ObjectIterator o_it(PluginFunctions::TARGET_OBJECTS, DATA_TETRAHEDRAL_MESH);
         o_it != PluginFunctions::objectsEnd(); ++o_it) {
        auto tet_obj = PluginFunctions::tetrahedralMeshObject(*o_it);
        auto tetmesh = tet_obj->mesh();

        tet_obj->materialNode()->set_point_size(12);

        tet_obj->meshNode()->drawMode(
                    ACG::SceneGraph::DrawModes::WIREFRAME
                  | ACG::SceneGraph::DrawModes::POINTS_COLORED
                  | ACG::SceneGraph::DrawModes::SOLID_FACES_COLORED|
                    DrawModes::Cells_flat() |
                    DrawModes::Vertices());

        tet_obj->materialNode()->enable_alpha_test(0.8);


        emit updatedObject(tet_obj->id(), UPDATE_COLOR);
    }
}

void MasterPlugin::slot_displace_constraint_vertex(){
    ACG::Vec3d displacement(tool_->displacementX->value(), tool_->displacementY->value(), tool_->displacementZ->value());
    // Tri meshes
    for (PluginFunctions::ObjectIterator o_it(PluginFunctions::TARGET_OBJECTS, DATA_TRIANGLE_MESH);
         o_it != PluginFunctions::objectsEnd(); ++o_it) {
        auto *tri_obj = PluginFunctions::triMeshObject(*o_it);
        auto *trimesh = tri_obj->mesh();

        tri_obj->materialNode()->set_point_size(12);

        tri_obj->meshNode()->drawMode(
                    ACG::SceneGraph::DrawModes::WIREFRAME
                  | ACG::SceneGraph::DrawModes::POINTS_COLORED
                  | ACG::SceneGraph::DrawModes::SOLID_FACES_COLORED);

        tri_obj->materialNode()->enable_alpha_test(0.8);


        VertexDisplacement::displace(mesh_, displacement, constraint_vhs_, false);

        TriangleLoop loop(mesh_, worldMesh_, q_min_, constraint_vhs_);

        loop.loop();

        std::cout << "Loop ended" << std::endl;

        std::cout << "Constraint vertices" << std::endl;
        for(auto v: constraint_vhs_){
            std::cout << v.first <<"," ;
        }
        std::cout << "\nTotal: " << constraint_vhs_.size() << std::endl;
        TriMesh newMesh = *trimesh;
        newMesh.clear();
        newMesh = mesh_;
        *trimesh = newMesh;
        trimesh->garbage_collection();

        if(tool_->showQualityCheckbox->isChecked())
            Highlight::highlight_triangles(*trimesh);
        emit updatedObject(tri_obj->id(), UPDATE_ALL);
    }
    // Tet meshes
    for (PluginFunctions::ObjectIterator o_it(PluginFunctions::TARGET_OBJECTS, DATA_TETRAHEDRAL_MESH);
         o_it != PluginFunctions::objectsEnd(); ++o_it) {
        auto *tet_obj = PluginFunctions::tetrahedralMeshObject(*o_it);
        auto *tetmesh = tet_obj->mesh();

        tet_obj->meshNode()->drawMode(
                    DrawModes::Cells_flat() |
                    DrawModes::Vertices());

        VertexDisplacement::displace(tetmesh_, displacement, constraint_vhs_);

        TetLoop loop3d(tetmesh_, q_min_, constraint_vhs_);
        loop3d.loop();

        std::cout << "Loop ended" << std::endl;

        *tetmesh = tetmesh_;
        tetmesh->collect_garbage();

        emit updatedObject(tet_obj->id(), UPDATE_ALL);
    }

}


void MasterPlugin::slot_start_experiment(){
    t_ = 0;
    timesteps_ = tool_->timestepsSpinBox->value();
    q_min_ = tool_->qualitySpinBox->value();
    bool useWorld = tool_->worldSpaceCheck->isChecked();
    bool autocomplete(true);

    tool_->nextButton->setEnabled(true);
    tool_->completeButton->setEnabled(true);
    tool_->beginExpButton->setEnabled(false);
    tool_->experiment_tabs->setEnabled(false);
    tool_->timestepsSpinBox->setEnabled(false);
    tool_->qualitySpinBox->setEnabled(false);
    tool_->experiment2D->setEnabled(false);
    tool_->experiment3D->setEnabled(false);
    tool_->splitCheckBox->setEnabled(false);
    tool_->splitSpinBox->setEnabled(false);
    tool_->worldSpaceCheck->setEnabled(false);
    tool_->angleSpinBox->setEnabled(false);

    // -------------------------- 2D --------------------
    for (PluginFunctions::ObjectIterator o_it(PluginFunctions::TARGET_OBJECTS, DATA_TRIANGLE_MESH);
         o_it != PluginFunctions::objectsEnd(); ++o_it) {
        auto *tri_obj = PluginFunctions::triMeshObject(*o_it);
        auto *trimesh = tri_obj->mesh();

        tri_obj->materialNode()->set_point_size(3.0);

        std::map<int,ACG::Vec3d> basePoints;

        for(auto vh: mesh_.vertices()){
            basePoints[vh.idx()] = mesh_.point(vh);
        }
        experiment2D_ = new Experiment(mesh_, worldMesh_,basePoints, q_min_, constraint_vhs_);
        std::cout << "Selected experiment: " << tool_->experiment2D->currentText().toStdString() << std::endl;

        slot_experiment_loop();
    }
    // -------------------------- 3D --------------------
    for (PluginFunctions::ObjectIterator o_it(PluginFunctions::TARGET_OBJECTS, DATA_TETRAHEDRAL_MESH);
         o_it != PluginFunctions::objectsEnd(); ++o_it) {
        auto *tet_obj = PluginFunctions::tetrahedralMeshObject(*o_it);
        auto *tetmesh = tet_obj->mesh();

        tet_obj->materialNode()->set_point_size(3.0);


        auto begin = std::chrono::high_resolution_clock::now();
        experiment3D_ = new Experiment3D(tetmesh_, -q_min_, constraint_vhs_, useWorld);
        std::cout << "Selected experiment: " <<
                     tool_->experiment3D->currentText().toStdString() << std::endl;
        if(autocomplete){
            slot_complete_experiment();
            TetLoop::displayIterationTime(begin, "Complete experiment time");
        }
    }
}

void MasterPlugin::slot_experiment_loop(){

    // -------------------------- 2D --------------------
    for (PluginFunctions::ObjectIterator o_it(PluginFunctions::TARGET_OBJECTS, DATA_TRIANGLE_MESH);
         o_it != PluginFunctions::objectsEnd(); ++o_it) {
        auto *tri_obj = PluginFunctions::triMeshObject(*o_it);
        auto *trimesh = tri_obj->mesh();

        /* indices 2D :
         * 0: Compress
         * 1: Stretch
        */
        switch (tool_->experiment2D->currentIndex()) {
            case 0:
                experiment2D_->compress2D(timesteps_, ++t_);
                break;
            case 1:
                experiment2D_->stretch2D(timesteps_, ++t_);
                break;
            default:
                break;
        }

        TriMesh newMesh = *trimesh;
        newMesh.clear();
        newMesh = mesh_;
        *trimesh = newMesh;
        trimesh->garbage_collection();

        if(tool_->showQualityCheckbox->isChecked())
            Highlight::highlight_triangles(*trimesh);

        emit updatedObject(tri_obj->id(), UPDATE_ALL);

        // reset when done with experiment
        if(t_ >= timesteps_){
            tool_->nextButton       ->setEnabled(false);
            tool_->completeButton   ->setEnabled(false);
            tool_->beginExpButton   ->setEnabled(true);
            tool_->experiment_tabs  ->setEnabled(true);
            tool_->timestepsSpinBox ->setEnabled(true);
            tool_->experiment2D     ->setEnabled(true);
            tool_->experiment3D     ->setEnabled(true);
            tool_->splitCheckBox    ->setEnabled(true);
            tool_->splitSpinBox     ->setEnabled(true);
            tool_->worldSpaceCheck  ->setEnabled(true);
            tool_->angleSpinBox     ->setEnabled(true);
            tool_->qualitySpinBox   ->setEnabled(true);
        }
    }
    // -------------------------- 3D --------------------
    for (PluginFunctions::ObjectIterator o_it(PluginFunctions::TARGET_OBJECTS, DATA_TETRAHEDRAL_MESH);
         o_it != PluginFunctions::objectsEnd(); ++o_it) {
        auto *tet_obj = PluginFunctions::tetrahedralMeshObject(*o_it);
        auto *tetmesh = tet_obj->mesh();
        double nbLoops = tool_->angleSpinBox->value();
        std::cout << "\033[1;35mTimestep: "<< t_<< "/"<< timesteps_<< "\033[0m" << std::endl;
        experiment3D_->generate_torsion_mesh(nbLoops/timesteps_);
        ++t_;

        *tetmesh = tetmesh_;
        tetmesh->collect_garbage();

        emit updatedObject(tet_obj->id(), UPDATE_ALL);

        // reset when done with experiment
        if(t_ >= timesteps_){
            double average(0.);
            TetLoop::Tet
                    worstInside(CellHandle(-1), 0),
                    worstBoundary(CellHandle(-1), 0);
            for(auto ch: tetmesh_.cells()){
                double quality = QualityEvaluation::evaluate(ch, tetmesh_);
                average += quality;
                if(tetmesh_.is_boundary(ch) && quality < worstBoundary.quality_){
                    worstBoundary.cell_handle_ = ch;
                    worstBoundary.quality_ = quality;
                }else if(!tetmesh_.is_boundary(ch) && quality < worstInside.quality_){
                    worstInside.cell_handle_ = ch;
                    worstInside.quality_ = quality;
                }
            }
            average /= tetmesh_.n_logical_cells();

            std::cout << "Experiment ended" << std::endl;
            std::cout << "Worst boundary cell: "<< worstBoundary.cell_handle_ <<
                      ", q="<< worstBoundary.quality_<< std::endl;
            std::cout << "Worst interior cell: "<< worstInside.cell_handle_ <<
                      ", q="<< worstInside.quality_<< std::endl;
            std::cout << "Average mesh quality: "<< average << std::endl;
            std::cout << "Experiment ended" << std::endl;
            tool_->nextButton       ->setEnabled(false);
            tool_->completeButton   ->setEnabled(false);
            tool_->beginExpButton   ->setEnabled(true);
            tool_->experiment_tabs  ->setEnabled(true);
            tool_->timestepsSpinBox ->setEnabled(true);
            tool_->experiment2D     ->setEnabled(true);
            tool_->experiment3D     ->setEnabled(true);
            tool_->splitCheckBox    ->setEnabled(true);
            tool_->splitSpinBox     ->setEnabled(true);
            tool_->worldSpaceCheck  ->setEnabled(true);
            tool_->angleSpinBox     ->setEnabled(true);
            tool_->qualitySpinBox   ->setEnabled(true);
        }
    }
}

// -------------------- Simple helpers ----------------------------
void MasterPlugin::highlight_constraints_vertices(TriMesh* _trimesh){
    for(auto vh: _trimesh->vertices()){
        if(!_trimesh->status(vh).deleted())
            _trimesh->set_color(vh, ACG::Vec4f(1,1,1,0));
    }
    for(auto vh_id: constraint_vhs_){
        auto vh = _trimesh->vertex_handle(vh_id.first);
        if(!_trimesh->status(vh).deleted())
            _trimesh->set_color(vh, ACG::Vec4f(1,0,0,1));
    }
}


// -------------------- Mesh generation ---------------------------
void MasterPlugin::slot_generate_base_mesh(){
    std::cout << "Generating a " << tool_->meshDimension->value() << "x" << tool_->meshDimension->value() << " "
              <<tool_->meshType->currentText().toStdString() << std::endl;
    switch (tool_->meshType->currentIndex()) {
    case 0:
        generate_tet_mesh();
        break;
    case 1:
        generate_triangular_mesh();
        break;
    default:
        generate_triangular_mesh();
        break;
    }

}
void MasterPlugin::generate_triangular_mesh(){

    int mesh_obj_id;
    emit addEmptyObject(DATA_TRIANGLE_MESH, mesh_obj_id);
    auto *mesh_obj = PluginFunctions::triMeshObject(mesh_obj_id);
    mesh_obj->setName("Mesh");
    mesh_obj->materialNode()->set_point_size(6.0);
    // Create vertices
    CustomMesh::VertexHandle vhandle[121];
    int count = 0;
    int dimension = tool_->meshDimension->value();

    auto mesh = mesh_obj->mesh();

    for (int i = 0; i < dimension+1; ++i) {
       for (int j = 0; j < dimension+1; ++j) {
           ACG::Vec3d p(0.2*i, 0.2*j, 0);
           vhandle[count] = mesh->add_vertex(p);
           mesh->set_color(vhandle[count++], ACG::Vec4f(1,1,1,0));
       }
    }

    // Create faces
    for (int i = 0; i < dimension; ++i) {
       for (int j = 0; j < dimension; ++j) {
           // First triangle
           CustomMesh::FaceHandle fh1 = mesh->add_face(
                        vhandle[i*(dimension+1)+j],
                        vhandle[i*(dimension+1)+j+1],
                        vhandle[(i+1)*(dimension+1)+j+1]);
           // Second triangle
           CustomMesh::FaceHandle fh2 = mesh->add_face(
                        vhandle[(i+1)*(dimension+1)+j+1],
                        vhandle[(i+1)*(dimension+1)+j],
                        vhandle[i*(dimension+1)+j]);
           mesh->set_color(fh1, ACG::Vec4f(1,1,1,1));
           mesh->set_color(fh2, ACG::Vec4f(1,1,1,1));
       }
    }

    mesh_obj->meshNode()->drawMode(
                ACG::SceneGraph::DrawModes::WIREFRAME
              | ACG::SceneGraph::DrawModes::POINTS_COLORED
              | ACG::SceneGraph::DrawModes::SOLID_FACES_COLORED);

    mesh_obj->materialNode()->enable_alpha_test(0.8);

    //  Request required status flags
    mesh->request_vertex_status();
    mesh->request_edge_status();
    mesh->request_face_status();
    emit updatedObject(mesh_obj->id(), UPDATE_ALL);

    mesh_ = *mesh;

}

void MasterPlugin::generate_tet_mesh(){
    int dimension = tool_->meshDimension->value();
    int mesh_obj_id;
    emit addEmptyObject(DATA_TETRAHEDRAL_MESH, mesh_obj_id);
    auto *mesh_obj = PluginFunctions::tetrahedralMeshObject(mesh_obj_id);
    mesh_obj->setName("Mesh");
    mesh_obj->materialNode()->set_point_size(6.0);

    // Create a mesh object
    auto mesh = mesh_obj->mesh();

    TetrahedralizedVoxelGridGenerator<TetrahedralMesh>::
            generate_mesh(dimension, dimension,dimension*2, *mesh);

    //split boundary edges
    if(tool_->splitCheckBox->isChecked()){
        int splits = tool_->splitSpinBox->value();
        std::cout<<" - running "<<splits<<" split passes on the boundary"<<std::endl;
        for(int i(0); i<splits; i++){
            //std::vector<EdgeHandle> to_split;
            EntityPriorityQueue<EdgeHandle, double> to_split;
            for(auto e: mesh->edges()){
                if(mesh->is_boundary(e)){
                    //to_split.push_back(e);
                    double len = (mesh->vertex(mesh->edge(e).from_vertex()) - mesh->vertex(mesh->edge(e).to_vertex())).norm();
                    to_split.push(e, -len);
                }
            }
            while(!to_split.empty()){
                auto e = to_split.front();
                to_split.pop();
                mesh->split_edge(e);
            }
            /*for(auto e: to_split){
                    mesh.split_edge(e);
                }*/
        }
    }

    QualityEvaluation::scaleMesh(*mesh);

    mesh_obj->meshNode()->drawMode(
                DrawModes::Cells_flat() |
                DrawModes::Vertices());

    mesh_obj->materialNode()->enable_alpha_test(0.8);
    emit updatedObject(mesh_obj->id(), UPDATE_ALL);

    tetmesh_ = *mesh;
}

TriMesh MasterPlugin::gen_world_mesh(){
    CustomMesh worldMesh;
    CustomMesh::VertexHandle vhandle[121];
    int count = 0;
    int dimension = tool_->meshDimension->value();
    for (int i = 0; i < dimension+1; ++i) {
       for (int j = 0; j < dimension+1; ++j) {
           ACG::Vec3d p(0.2*i, 0.2*j, 0);
           vhandle[count] = worldMesh.add_vertex(p);
       }
    }
    // Create faces
    for (int i = 0; i < dimension; ++i) {
       for (int j = 0; j < dimension; ++j) {
           // First triangle
           CustomMesh::FaceHandle fh1 = worldMesh.add_face(
                        worldMesh.vertex_handle(i*(dimension+1)+j),
                        worldMesh.vertex_handle(i*(dimension+1)+j+1),
                        worldMesh.vertex_handle((i+1)*(dimension+1)+j+1));
           // Second triangle
           CustomMesh::FaceHandle fh2 = worldMesh.add_face(
                        worldMesh.vertex_handle((i+1)*(dimension+1)+j+1),
                        worldMesh.vertex_handle((i+1)*(dimension+1)+j),
                        worldMesh.vertex_handle(i*(dimension+1)+j));
       }
    }
    //  Request required status flags
    worldMesh.request_vertex_status();
    worldMesh.request_edge_status();
    worldMesh.request_face_status();
    return worldMesh;
}


#if QT_VERSION < 0x050000
Q_EXPORT_PLUGIN2(masterplugin, MasterPlugin);
#endif
