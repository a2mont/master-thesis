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

#include <ObjectTypes/TriangleMesh/TriangleMesh.hh>
#include <OpenVolumeMesh/Geometry/VectorT.hh>
#include <OpenVolumeMesh/Mesh/PolyhedralMesh.hh>
#include <ObjectTypes/TetrahedralMesh/TetrahedralMesh.hh>

#include "MainLoop.hh"
#include "MasterToolbar.hh"
#include "Highlight.hh"

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
    void log(Logtype _type, QString _message);
    void log(QString _message);

    // LoadSaveInterface
    void addEmptyObject(DataType _type, int& _id);
    void updatedObject(int _identifier, const UpdateType& _type);

    // ToolboxInterface
    void addToolbox(QString _name, QWidget* _widget);
    //PickingInterface
    void addPickMode(const std::string& _mode);
    void addHiddenPickMode(const std::string& _mode);


private slots:
    // initialization functions
    void initializePlugin();
    void pluginsInitialized();

    void slotMouseEvent(QMouseEvent* _event) ;
    void slot_show_constraint_vertex();
    void slot_displace_constraint_vertex();
    void slot_generate_base_mesh();
    void slot_show_quality();
    void slot_clear_constraints();


public :
  
  ~MasterPlugin() {};
  
  QString name() { return QString("Master Plugin"); };
  
  QString description() { return QString("Plugin used for my thesis !"); };
  

private:
    typedef OpenMesh::TriMesh_ArrayKernelT<> CustomMesh;
    typedef OpenVolumeMesh::GeometricPolyhedralMeshV3f MyMesh;
    // The toolbox widget and the button in it
    MasterToolbar* tool_;

    ACG::HaltonColors hcolors_;

    //store selected vertices
    std::map<int,int> constraint_vhs_;

    const double q_min_ = 0.5;

    TriMesh trimesh_;
    TetrahedralMesh tetmesh_;
    void generate_triangular_mesh();
    void generate_tet_mesh();
};
#endif
