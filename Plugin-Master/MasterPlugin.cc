#include "MasterPlugin.hh"

#include <ObjectTypes/PolyMesh/PolyMesh.hh>
#include <ObjectTypes/TriangleMesh/TriangleMesh.hh>
#include "OpenFlipper/BasePlugin/PluginFunctions.hh"



void MasterPlugin::initializePlugin()
{  
    tool_ = new QWidget();
    QSize size(600, 300);
    tool_->resize(size);

    // Mesh Generation
    dimensionComboBox_ = new QComboBox();
    dimensionComboBox_->addItem("Triangular Mesh");
    dimensionComboBox_->addItem("Tetrahedral Mesh");
    auto generateButton = new QPushButton("Generate Mesh");

    // Vertex selection
    QLabel* labelVertex = new QLabel("Pick constraint vertex");
    pickButton_ = new QPushButton("Select vertex");
    pickButton_->setCheckable(true);
    auto generationLabel = new QLabel("Generate base mesh");
    showQualityCheckBox_ = new QCheckBox("Show quality");
    showQualityCheckBox_->setChecked(true);

    // Vertext displacement
    QLabel* labelDisplacement = new QLabel("Displacement values for constrained vertex");
    xValue_ = new QDoubleSpinBox();
    xValue_->setPrefix("x: ");
    xValue_->setSingleStep(0.01);
    xValue_->setRange(-1.0, 1.0);
    xValue_->setValue(.05);
    yValue_ = new QDoubleSpinBox();
    yValue_->setPrefix("y: ");
    yValue_->setValue(0.0);
    yValue_->setSingleStep(0.01);
    yValue_->setRange(-1.0, 1.0);
    zValue_ = new QDoubleSpinBox();
    zValue_->setPrefix("z: ");
    zValue_->setValue(0.0);
    zValue_->setSingleStep(0.01);
    zValue_->setRange(-1.0, 1.0);

    displaceButton_ = new QPushButton("Displace vertex");


    // Layout
    QGridLayout* grid = new QGridLayout();
    grid->addWidget(generationLabel, 0,0);
    grid->addWidget(dimensionComboBox_, 2,0);
    grid->addWidget(generateButton, 2,1);
    grid->addWidget(labelVertex, 3, 0);
    grid->addWidget(pickButton_, 4, 0);
    grid->addWidget(showQualityCheckBox_, 4,1);
    grid->addWidget(labelDisplacement, 5,0);
    grid->addWidget(xValue_, 6,0);
    grid->addWidget(yValue_, 6,1);
    grid->addWidget(zValue_, 6,2);
    grid->addWidget(displaceButton_, 7,0);
    tool_->setLayout(grid);

    // Connections
    connect(generateButton, SIGNAL(clicked()), this, SLOT(slot_generate_base_mesh()));
    connect(pickButton_, SIGNAL(clicked()), this, SLOT(slot_show_constraint_vertex()));
    connect(displaceButton_,  SIGNAL(clicked()), this, SLOT(slot_displace_constraint_vertex()));
    connect(showQualityCheckBox_, SIGNAL(toggled(bool)), this, SLOT(slot_show_quality()));

    // Arbitrary id for constraint vertex
    constraint_vh_ = 15;


    // Add the Toolbox
    emit addToolbox("Master", tool_);
}

void MasterPlugin::pluginInitialized(){
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
                        //set constraint vertices (for #selected vertices >1)
//                        for(int i=0; i<2; ++i)
//                            if(tool_->vertex_number_cb->currentIndex() == i)
                        constraint_vh_ = vh.idx();

                        slot_show_constraint_vertex();

                        return;
                    }
                }
            }
        }
    }

    emit updateView();
}
void MasterPlugin::slot_show_quality(){
    for (PluginFunctions::ObjectIterator o_it(PluginFunctions::TARGET_OBJECTS, DATA_TRIANGLE_MESH);
         o_it != PluginFunctions::objectsEnd(); ++o_it) {
        auto tri_obj = PluginFunctions::triMeshObject(*o_it);
        auto trimesh = tri_obj->mesh();


        if(showQualityCheckBox_->isChecked()){
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
    if(pickButton_->isChecked()) {
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

        trimesh->set_color(trimesh->vertex_handle(constraint_vh_), ACG::Vec4f(1,0,0,1));

        tri_obj->meshNode()->drawMode(
                    ACG::SceneGraph::DrawModes::WIREFRAME
                  | ACG::SceneGraph::DrawModes::POINTS_COLORED
                  | ACG::SceneGraph::DrawModes::SOLID_FACES_COLORED);

        tri_obj->materialNode()->enable_alpha_test(0.8);


        emit updatedObject(tri_obj->id(), UPDATE_COLOR);
    }
}

void MasterPlugin::slot_displace_constraint_vertex(){
    ACG::Vec3d displacement(xValue_->value(), yValue_->value(), zValue_->value());
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


        MainLoop loop(trimesh_, q_min_);


        loop.loop(displacement, constraint_vh_, false);
        std::cout << "Loop ended" << std::endl;
        *trimesh = trimesh_;
        trimesh->garbage_collection();

        trimesh->set_color(trimesh->vertex_handle(constraint_vh_), ACG::Vec4f(1,0,0,1));

        if(showQualityCheckBox_->isChecked())
            Highlight::highlight_triangles(*trimesh);

        emit updatedObject(tri_obj->id(), UPDATE_ALL);
    }

}

void MasterPlugin::slot_generate_base_mesh(){
    std::cout << "Generating a " << dimensionComboBox_->currentText().toStdString() << std::endl;
    switch (dimensionComboBox_->currentIndex()) {
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

    for (int i = 0; i < 11; ++i) {
       for (int j = 0; j < 11; ++j) {
           ACG::Vec3d p(-1 + 0.2*i, -1 + 0.2*j, 0);
           vhandle[count] = mesh->add_vertex(p);
           mesh->set_color(vhandle[count++], ACG::Vec4f(1,1,1,0));
       }
    }

    // Create faces
    for (int i = 0; i < 10; ++i) {
       for (int j = 0; j < 10; ++j) {
           // First triangle
           CustomMesh::FaceHandle fh1 = mesh->add_face(vhandle[i*11+j], vhandle[i*11+j+1], vhandle[(i+1)*11+j+1]);
           // Second triangle
           CustomMesh::FaceHandle fh2 = mesh->add_face(vhandle[(i+1)*11+j+1], vhandle[(i+1)*11+j], vhandle[i*11+j]);
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
