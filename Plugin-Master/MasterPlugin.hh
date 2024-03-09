#ifndef MASTERPLUGIN_HH
#define MASTERPLUGIN_HH

#include <QObject>

#include <OpenFlipper/common/Types.hh>
#include <OpenFlipper/BasePlugin/BaseInterface.hh>
#include <OpenFlipper/BasePlugin/ToolboxInterface.hh>
#include <OpenFlipper/BasePlugin/LoggingInterface.hh>
#include <OpenFlipper/BasePlugin/LoadSaveInterface.hh>
#include <OpenFlipper/BasePlugin/MouseInterface.hh>
#include <OpenFlipper/BasePlugin/PickingInterface.hh>
#include <ACG/Utils/HaltonColors.hh>
#include <ACG/Scenegraph/LineNode.hh>

#include <unistd.h>
#include <ObjectTypes/TriangleMesh/TriangleMesh.hh>
#include <OpenVolumeMesh/Geometry/VectorT.hh>
#include <OpenVolumeMesh/Mesh/PolyhedralMesh.hh>
#include <ObjectTypes/TetrahedralMesh/TetrahedralMesh.hh>

#include "TriangleLoop.hh"
#include "TetLoop.hh"
#include "MasterToolbar.hh"
#include "DrawModes.hh"
#include "Highlight.hh"
#include "Experiments.hh"
#include "Experiments3D.hh"
#include "TetrahedralizedVoxelGridGenerator.hh"
#include "Tests.hh"

class MasterPlugin : public QObject, BaseInterface, ToolboxInterface, LoggingInterface, LoadSaveInterface, MouseInterface, PickingInterface
{
Q_OBJECT
    Q_INTERFACES(BaseInterface)
    Q_INTERFACES(ToolboxInterface)
    Q_INTERFACES(LoggingInterface)
    Q_INTERFACES(LoadSaveInterface)
    Q_INTERFACES(MouseInterface)
    Q_INTERFACES(PickingInterface)

#if QT_VERSION >= 0x050000
    Q_PLUGIN_METADATA(IID "org.OpenFlipper.Plugins.Plugin-Master")
#endif
signals:
    void updateView();


    //LoggingInterface
    void log(Logtype _type, QString _message) override;
    void log(QString _message) override;

    // LoadSaveInterface
    void addEmptyObject(DataType _type, int& _id) override;
    void updatedObject(int _identifier, const UpdateType& _type) override;

    // ToolboxInterface
    void addToolbox(QString _name, QWidget* _widget, QIcon* _icon) override;
    //PickingInterface
    void addPickMode(const std::string& _mode) override;
    void addHiddenPickMode(const std::string& _mode) override;


private slots:
    // initialization functions
    void initializePlugin() override;
    void pluginsInitialized() override;

    void slotMouseEvent(QMouseEvent* _event)  override;
    void slot_show_constraint_vertex();
    void slot_displace_constraint_vertex();
    void slot_generate_base_mesh();
    void slot_show_quality();
    void slot_clear_constraints();
    void slot_start_experiment();
    void slot_experiment_loop();
    void slot_complete_experiment();


public :

  ~MasterPlugin() {};

  QString name() override { return QString("Master Plugin"); };

  QString description() override { return QString("Plugin used for my thesis !"); };

private:
    typedef OpenMesh::TriMesh_ArrayKernelT<OpenMesh::DefaultTraitsDouble> CustomMesh;
    typedef OpenVolumeMesh::GeometricPolyhedralMeshV3f MyMesh;
    // The toolbox widget and the button in it
    MasterToolbar* tool_;

    ACG::HaltonColors hcolors_;

    //store selected vertices
    std::map<int,int> constraint_vhs_;

    double q_min_;

    int timesteps_ = 0;
    int t_ = 0;
    TriMesh mesh_;
    TriMesh worldMesh_;
    TetrahedralMesh tetmesh_;
    Experiment *experiment2D_;
    Experiment3D* experiment3D_;

    void generate_triangular_mesh();
    void generate_tet_mesh();
    void highlight_constraints_vertices(TriMesh *_trimesh);
    TriMesh gen_world_mesh();
};
#endif
