#include "MasterPlugin.hh"

#include <ObjectTypes/PolyMesh/PolyMesh.hh>
#include <ObjectTypes/TriangleMesh/TriangleMesh.hh>
#include "OpenFlipper/BasePlugin/PluginFunctions.hh"



void MasterPlugin::initializePlugin()
{  
    tool_ = new QWidget();
    QSize size(600, 300);
    tool_->resize(size);

    auto generateButton = new QPushButton(tr("Generate Mesh"));
    showQualityCheckBox_ = new QCheckBox(tr("Show quality"));
    showQualityCheckBox_->setChecked(true);

    pickButton_ = new QPushButton(tr("Select vertex"));
    pickButton_->setCheckable(true);
    QLabel* labelVertex = new QLabel("Pick constraint vertex");


    QLabel* labelDisplacement = new QLabel("Displacement values for constrained vertex");
    xValue_ = new QDoubleSpinBox();
    xValue_->setPrefix(tr("x: "));
    xValue_->setSingleStep(0.01);
    xValue_->setRange(-1.0, 1.0);
    xValue_->setValue(.05);
    yValue_ = new QDoubleSpinBox();
    yValue_->setPrefix(tr("y: "));
    yValue_->setValue(.05);
    yValue_->setSingleStep(0.01);
    yValue_->setRange(-1.0, 1.0);
    zValue_ = new QDoubleSpinBox();
    zValue_->setPrefix(tr("z: "));
    zValue_->setValue(0.0);
    zValue_->setSingleStep(0.01);
    zValue_->setRange(-1.0, 1.0);

    displaceButton_ = new QPushButton(tr("Displace vertex"));

    QGridLayout* grid = new QGridLayout();
    grid->addWidget(generateButton, 0,0);
    grid->addWidget(showQualityCheckBox_, 0,1);
    grid->addWidget(labelVertex, 1, 0);
    grid->addWidget(pickButton_, 2, 0);
    grid->addWidget(labelDisplacement, 3,0);
    grid->addWidget(xValue_, 4,0);
    grid->addWidget(yValue_, 4,1);
    grid->addWidget(zValue_, 4,2);
    grid->addWidget(displaceButton_, 5,0);
    tool_->setLayout(grid);

    // Connect button to slotButtonClicked()
    connect(generateButton, SIGNAL(clicked()), this, SLOT(slot_generate_base_mesh()));
    connect(pickButton_, SIGNAL(clicked()), this, SLOT(slot_show_constraint_vertex()));
    connect(displaceButton_,  SIGNAL(clicked()), this, SLOT(slot_displace_constraint_vertex()));
    connect(showQualityCheckBox_, SIGNAL(toggled(bool)), this, SLOT(slot_show_quality()));

    constraint_vh_ = 25;


    // Add the Toolbox
    emit addToolbox(tr("Master"), tool_);
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

        for(auto vh : trimesh->vertices()) {
            if(constraint_vh_ == vh.idx())
                trimesh->set_color(vh, ACG::Vec4f(1,0,0,1));
            else
                trimesh->set_color(vh, ACG::Vec4f(1,1,1,0));
        }

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

        for(auto fh: trimesh->faces()){
            trimesh->set_color(fh, ACG::Vec4f(1,1,1,1));
        }
        tri_obj->meshNode()->drawMode(
                    ACG::SceneGraph::DrawModes::WIREFRAME
                  | ACG::SceneGraph::DrawModes::POINTS_COLORED
                  | ACG::SceneGraph::DrawModes::SOLID_FACES_COLORED);

        tri_obj->materialNode()->enable_alpha_test(0.8);


        MainLoop loop(mesh_, q_min_);


        loop.loop(displacement, constraint_vh_, false);
        std::cout << "Loop ended" << std::endl;
        *trimesh = mesh_;
        trimesh->garbage_collection();
        if(showQualityCheckBox_->isChecked())
            Highlight::highlight_triangles(*trimesh);

        emit updatedObject(tri_obj->id(), UPDATE_ALL);
    }
}

void MasterPlugin::slot_generate_base_mesh(){
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

    mesh_ = *mesh;

}


#if QT_VERSION < 0x050000
Q_EXPORT_PLUGIN2(masterplugin, MasterPlugin);
#endif
