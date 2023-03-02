#include "MasterPlugin.hh"

#include <ObjectTypes/PolyMesh/PolyMesh.hh>
#include <ObjectTypes/TriangleMesh/TriangleMesh.hh>
#include "OpenFlipper/BasePlugin/PluginFunctions.hh"



void MasterPlugin::initializePlugin()
{  
    tool_ = new QWidget();
    QSize size(300, 300);
    tool_->resize(size);

    // Create button that can be toggled
    // to (de)activate plugin's picking mode
    pickButton_ = new QPushButton(tr("Select object"));
    pickButton_->setCheckable(true);
    QLabel* label = new QLabel("Pick constraint vertex");

    QGridLayout* grid = new QGridLayout();
    grid->addWidget(label, 0, 0);
    grid->addWidget(pickButton_, 1, 0);
    tool_->setLayout(grid);

    // Connect button to slotButtonClicked()
    connect(pickButton_, SIGNAL(clicked()), this, SLOT(slot_show_constraint_vertex()));
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
                        //set constraint vertices
//                        for(int i=0; i<2; ++i)
//                            if(tool_->vertex_number_cb->currentIndex() == i)
                        constraint_vhs_ = vh.idx();

                        slot_show_constraint_vertex();

                        return;
                    }
                }
            }
        }
    }

    emit updateView();
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
            if(constraint_vhs_ == vh.idx())
                trimesh->set_color(vh, ACG::Vec4f(1,0,0,1));
            else
                trimesh->set_color(vh, ACG::Vec4f(1,1,1,0));
        }

        tri_obj->meshNode()->drawMode(ACG::SceneGraph::DrawModes::WIREFRAME
        | ACG::SceneGraph::DrawModes::SOLID_SMOOTH_SHADED | ACG::SceneGraph::DrawModes::POINTS_COLORED);

        tri_obj->materialNode()->enable_alpha_test(0.8);

        emit updatedObject(tri_obj->id(), UPDATE_COLOR);
    }
}


#if QT_VERSION < 0x050000
Q_EXPORT_PLUGIN2(masterplugin, MasterPlugin);
#endif
