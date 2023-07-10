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


    // Add the Toolbox
    emit addToolbox("Master", tool_);
}

void MasterPlugin::pluginsInitialized(){
    // Arbitrary id for constraint vertex
    constraint_vhs_[51] = 51;
    slot_generate_base_mesh();
    slot_show_quality();
    slot_show_constraint_vertex();
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
        for(auto vh: trimesh->vertices())
            trimesh->set_color(vh, ACG::Vec4f(1,1,1,0));
        for(auto vh: constraint_vhs_)
            trimesh->set_color(trimesh->vertex_handle(vh.first), ACG::Vec4f(1,0,0,1));

        tri_obj->meshNode()->drawMode(
                    ACG::SceneGraph::DrawModes::WIREFRAME
                  | ACG::SceneGraph::DrawModes::POINTS_COLORED
                  | ACG::SceneGraph::DrawModes::SOLID_FACES_COLORED);

        tri_obj->materialNode()->enable_alpha_test(0.8);


        emit updatedObject(tri_obj->id(), UPDATE_COLOR);
    }
}

void MasterPlugin::slot_displace_constraint_vertex(){
    ACG::Vec3d displacement(tool_->displacementX->value(), tool_->displacementY->value(), tool_->displacementZ->value());
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


        MainLoop loop(trimesh_, q_min_, constraint_vhs_);

        loop.loop(displacement, false);
        std::cout << "Loop ended" << std::endl;
        *trimesh = trimesh_;
        trimesh->garbage_collection();

        for(auto vh: constraint_vhs_)
            trimesh->set_color(trimesh->vertex_handle(vh.first), ACG::Vec4f(1,0,0,1));

        if(tool_->showQualityCheckbox->isChecked())
            Highlight::highlight_triangles(*trimesh);

        emit updatedObject(tri_obj->id(), UPDATE_ALL);

    }

}

void MasterPlugin::slot_generate_base_mesh(){
    std::cout << "Generating a " << tool_->meshDimension->value() << "x" << tool_->meshDimension->value() << " "
              <<tool_->meshType->currentText().toStdString() << std::endl;
    switch (tool_->meshType->currentIndex()) {
    case 0:
        generate_triangular_mesh();
        break;
    case 1:
        generate_tet_mesh();
        break;
    default:
        generate_triangular_mesh();
        break;
    }

}

void MasterPlugin::slot_start_experiment(){
    for (PluginFunctions::ObjectIterator o_it(PluginFunctions::TARGET_OBJECTS, DATA_TRIANGLE_MESH);
         o_it != PluginFunctions::objectsEnd(); ++o_it) {
        auto *tri_obj = PluginFunctions::triMeshObject(*o_it);
        auto *trimesh = tri_obj->mesh();

        int timesteps = tool_->timestepsSpinBox->value();
        double pause = tool_->pauseSpinBox->value();
        std::cout << "Selected experiment: " << tool_->experiment2D->currentText().toStdString() << std::endl;
        if(tool_->experiment_tabs->currentIndex() == 0){
            for (int i = 0; i < timesteps; ++i) {
                /* indices 2D :
                 * 0: Bend
                 * 1: Compress
                 * 2: Stretch
                */
                switch (tool_->experiment2D->currentIndex()) {
                    case 0:
                        Experiments::bend2D(trimesh_, constraint_vhs_, timesteps, pause);
                        break;
                    case 1:
                        Experiments::compress2D(trimesh_, constraint_vhs_, timesteps, pause);
                        break;
                    case 2:
                        Experiments::stretch2D(trimesh_, constraint_vhs_, timesteps, pause);
                        break;
                    default:
                        break;
                }

            }

        }
        std::cout << "Experiment ended" << std::endl;
        *trimesh = trimesh_;
        trimesh->garbage_collection();

        if(tool_->showQualityCheckbox->isChecked())
            Highlight::highlight_triangles(*trimesh);
        emit updatedObject(tri_obj->id(), UPDATE_ALL);

    }
}

// -------------------- Mesh generation ---------------------------
void MasterPlugin::generate_triangular_mesh(){

    int mesh_obj_id;
    emit addEmptyObject(DATA_TRIANGLE_MESH, mesh_obj_id);
    auto *mesh_obj = PluginFunctions::triMeshObject(mesh_obj_id);
    mesh_obj->setName("Mesh");
    mesh_obj->materialNode()->set_point_size(6.0);
    // Create vertices
    CustomMesh::VertexHandle vhandle[121];
    int count = 0;

    auto mesh = mesh_obj->mesh();

    for (int i = 0; i < tool_->meshDimension->value()+1; ++i) {
       for (int j = 0; j < tool_->meshDimension->value()+1; ++j) {
           ACG::Vec3d p(-1 + 0.2*i, -1 + 0.2*j, 0);
           vhandle[count] = mesh->add_vertex(p);
           mesh->set_color(vhandle[count++], ACG::Vec4f(1,1,1,0));
       }
    }

    // Create faces
    for (int i = 0; i < tool_->meshDimension->value(); ++i) {
       for (int j = 0; j < tool_->meshDimension->value(); ++j) {
           // First triangle
           CustomMesh::FaceHandle fh1 = mesh->add_face(vhandle[i*(tool_->meshDimension->value()+1)+j],
                   vhandle[i*(tool_->meshDimension->value()+1)+j+1], vhandle[(i+1)*(tool_->meshDimension->value()+1)+j+1]);
           // Second triangle
           CustomMesh::FaceHandle fh2 = mesh->add_face(vhandle[(i+1)*(tool_->meshDimension->value()+1)+j+1],
                   vhandle[(i+1)*(tool_->meshDimension->value()+1)+j], vhandle[i*(tool_->meshDimension->value()+1)+j]);
           mesh->set_color(fh1, ACG::Vec4f(1,1,1,1));
           mesh->set_color(fh2, ACG::Vec4f(1,1,1,1));
       }
    }

    mesh_obj->meshNode()->drawMode(
                ACG::SceneGraph::DrawModes::WIREFRAME
              | ACG::SceneGraph::DrawModes::POINTS_COLORED
              | ACG::SceneGraph::DrawModes::SOLID_FACES_COLORED);

    mesh_obj->materialNode()->enable_alpha_test(0.8);
    emit updatedObject(mesh_obj->id(), UPDATE_ALL);

    trimesh_ = *mesh;
}


void MasterPlugin::generate_tet_mesh(){

    int mesh_obj_id;
    emit addEmptyObject(DATA_TETRAHEDRAL_MESH, mesh_obj_id);
    auto *mesh_obj = PluginFunctions::tetrahedralMeshObject(mesh_obj_id);
    mesh_obj->setName("Mesh");
    mesh_obj->materialNode()->set_point_size(6.0);

    // Create a mesh object
    auto mesh = mesh_obj->mesh();

    OpenVolumeMesh::VertexHandle v0 = mesh->add_vertex(ACG::Vec3d(-1.0, 0.0, 0.0));
    OpenVolumeMesh::VertexHandle v1 = mesh->add_vertex(ACG::Vec3d( 0.0, 0.0, 1.0));
    OpenVolumeMesh::VertexHandle v2 = mesh->add_vertex(ACG::Vec3d( 1.0, 0.0, 0.0));
    OpenVolumeMesh::VertexHandle v3 = mesh->add_vertex(ACG::Vec3d( 0.0, 0.0,-1.0));
    OpenVolumeMesh::VertexHandle v4 = mesh->add_vertex(ACG::Vec3d( 0.0, 1.0, 0.0));
    std::vector<OpenVolumeMesh::VertexHandle> vertices;
    // Add faces
    vertices.push_back(v0); vertices.push_back(v1);vertices.push_back(v4);
    OpenVolumeMesh::FaceHandle f0 = mesh->add_face(vertices);
    vertices.clear();
    vertices.push_back(v1); vertices.push_back(v2);vertices.push_back(v4);
    OpenVolumeMesh::FaceHandle f1 = mesh->add_face(vertices);
    vertices.clear();
    vertices.push_back(v0); vertices.push_back(v1);vertices.push_back(v2);
    OpenVolumeMesh::FaceHandle f2 = mesh->add_face(vertices);
    vertices.clear();
    vertices.push_back(v0); vertices.push_back(v4);vertices.push_back(v2);
    OpenVolumeMesh::FaceHandle f3 = mesh->add_face(vertices);
    vertices.clear();
    vertices.push_back(v0); vertices.push_back(v4);vertices.push_back(v3);
    OpenVolumeMesh::FaceHandle f4 = mesh->add_face(vertices);
    vertices.clear();
    vertices.push_back(v2); vertices.push_back(v3);vertices.push_back(v4);
    OpenVolumeMesh::FaceHandle f5 = mesh->add_face(vertices);
    vertices.clear();
    vertices.push_back(v0); vertices.push_back(v2);vertices.push_back(v3);
    OpenVolumeMesh::FaceHandle f6 = mesh->add_face(vertices);
    std::vector<OpenVolumeMesh::HalfFaceHandle> halffaces;
    // Add first tetrahedron
    halffaces.push_back(mesh->halfface_handle(f0, 1));
    halffaces.push_back(mesh->halfface_handle(f1, 1));
    halffaces.push_back(mesh->halfface_handle(f2, 0));
    halffaces.push_back(mesh->halfface_handle(f3, 1));
    mesh->add_cell(halffaces);
    // Add second tetrahedron
    halffaces.clear();
    halffaces.push_back(mesh->halfface_handle(f4, 1));
    halffaces.push_back(mesh->halfface_handle(f5, 1));
    halffaces.push_back(mesh->halfface_handle(f3, 0));
    halffaces.push_back(mesh->halfface_handle(f6, 0));
    mesh->add_cell(halffaces);


    mesh_obj->meshNode()->drawMode(
              ACG::SceneGraph::DrawModes::WIREFRAME
            | ACG::SceneGraph::DrawModes::POINTS_COLORED
            | ACG::SceneGraph::DrawModes::SOLID_FACES_COLORED);

    mesh_obj->materialNode()->enable_alpha_test(0.8);
    emit updatedObject(mesh_obj->id(), UPDATE_ALL);

    tetmesh_ = *mesh;

  return;


}

#if QT_VERSION < 0x050000
Q_EXPORT_PLUGIN2(masterplugin, MasterPlugin);
#endif
