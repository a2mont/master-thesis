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
    meshSize_ = new QSpinBox();
    meshSize_->setPrefix("Mesh size: ");
    meshSize_->setMinimum(1);
    meshSize_->setValue(10);
    auto generateButton = new QPushButton("Generate Mesh");

    // Vertex selection
    QLabel* labelVertex = new QLabel("Pick constraint vertices");
    pickButton_ = new QPushButton("Select mode");
    pickButton_->setCheckable(true);
    auto generationLabel = new QLabel("Generate base mesh");
    auto clearButton = new QPushButton("Clear vertices");
    showQualityCheckBox_ = new QCheckBox("Show quality");
    showQualityCheckBox_->setChecked(true);

    // Vertext displacement
    QLabel* labelDisplacement = new QLabel("Displacement values for\nconstrained vertices");
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

    displaceButton_ = new QPushButton("Displace vertices");


    // Layout
    QGridLayout* grid = new QGridLayout();
    grid->addWidget(generationLabel, 0,0);
    grid->addWidget(dimensionComboBox_, 1,0);
    grid->addWidget(meshSize_, 2,0);
    grid->addWidget(generateButton, 2,2);
    grid->addWidget(labelVertex, 3, 0);
    grid->addWidget(showQualityCheckBox_, 4,0);
    grid->addWidget(pickButton_, 5, 0);
    grid->addWidget(clearButton,5,2);
    grid->addWidget(labelDisplacement, 6,0);
    grid->addWidget(xValue_, 7,0);
    grid->addWidget(yValue_, 7,1);
    grid->addWidget(zValue_, 7,2);
    grid->addWidget(displaceButton_, 8,2);
    tool_->setLayout(grid);

    // Connections
    connect(generateButton, SIGNAL(clicked()), this, SLOT(slot_generate_base_mesh()));
    connect(pickButton_, SIGNAL(clicked()), this, SLOT(slot_show_constraint_vertex()));
    connect(displaceButton_,  SIGNAL(clicked()), this, SLOT(slot_displace_constraint_vertex()));
    connect(showQualityCheckBox_, SIGNAL(toggled(bool)), this, SLOT(slot_show_quality()));
    connect(clearButton, SIGNAL(clicked()), this, SLOT(slot_clear_constraints()));

    // Arbitrary id for constraint vertex
    constraint_vhs_[51] = 51;

    // Add the Toolbox
    emit addToolbox("Master", tool_);
}

void MasterPlugin::pluginsInitialized(){
    slot_generate_base_mesh();
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


        MainLoop loop(trimesh_, q_min_, constraint_vhs_);

        loop.loop(displacement, false);
        std::cout << "Loop ended" << std::endl;
        *trimesh = trimesh_;
        trimesh->garbage_collection();

        for(auto vh: constraint_vhs_)
            trimesh->set_color(trimesh->vertex_handle(vh.first), ACG::Vec4f(1,0,0,1));

        if(showQualityCheckBox_->isChecked())
            Highlight::highlight_triangles(*trimesh);

        emit updatedObject(tri_obj->id(), UPDATE_ALL);

    }

}

void MasterPlugin::slot_generate_base_mesh(){
    std::cout << "Generating a " << meshSize_->value() << "x" << meshSize_->value() << " "
              <<dimensionComboBox_->currentText().toStdString() << std::endl;
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

    for (int i = 0; i < meshSize_->value()+1; ++i) {
       for (int j = 0; j < meshSize_->value()+1; ++j) {
           ACG::Vec3d p(-1 + 0.2*i, -1 + 0.2*j, 0);
           vhandle[count] = mesh->add_vertex(p);
           mesh->set_color(vhandle[count++], ACG::Vec4f(1,1,1,0));
       }
    }

    // Create faces
    for (int i = 0; i < meshSize_->value(); ++i) {
       for (int j = 0; j < meshSize_->value(); ++j) {
           // First triangle
           CustomMesh::FaceHandle fh1 = mesh->add_face(vhandle[i*(meshSize_->value()+1)+j],
                   vhandle[i*(meshSize_->value()+1)+j+1], vhandle[(i+1)*(meshSize_->value()+1)+j+1]);
           // Second triangle
           CustomMesh::FaceHandle fh2 = mesh->add_face(vhandle[(i+1)*(meshSize_->value()+1)+j+1],
                   vhandle[(i+1)*(meshSize_->value()+1)+j], vhandle[i*(meshSize_->value()+1)+j]);
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
