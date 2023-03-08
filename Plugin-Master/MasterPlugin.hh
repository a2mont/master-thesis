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

#include "MasterToolbar.hh"
class MasterPlugin : public QObject, BaseInterface, ToolboxInterface, LoggingInterface, MouseInterface, PickingInterface
{
Q_OBJECT
    Q_INTERFACES(BaseInterface)
    Q_INTERFACES(ToolboxInterface)
    Q_INTERFACES(LoggingInterface)
    Q_INTERFACES(MouseInterface)
    Q_INTERFACES(PickingInterface)

#if QT_VERSION >= 0x050000
    Q_PLUGIN_METADATA(IID "org.OpenFlipper.Plugins.Plugin-Master")
#endif
signals:
    void updateView();
    void updatedObject(int _identifier, const UpdateType& _type);


    //LoggingInterface
    void log(Logtype _type, QString _message);
    void log(QString _message);

    // ToolboxInterface
    void addToolbox(QString _name, QWidget* _widget);
    //PickingInterface
    void addPickMode(const std::string& _mode);
    void addHiddenPickMode(const std::string& _mode);


private slots:
    // initialization functions
    void initializePlugin();
    void pluginInitialized();

    void slotMouseEvent(QMouseEvent* _event) ;
    void slot_show_constraint_vertex();
    void slot_displace_constraint_vertex();


public :
  
  ~MasterPlugin() {};
  
  QString name() { return QString("Master Plugin"); };
  
  QString description() { return QString("Plugin used for my thesis !"); };
  

private:
    // The toolbox widget and the button in it
    QWidget* tool_;
    QPushButton* pickButton_;
    QPushButton* displaceButton_;
    QDoubleSpinBox* xValue_;
    QDoubleSpinBox* yValue_;
    QDoubleSpinBox* zValue_;

    QSpinBox* iterationsSpinbox_;
    ACG::HaltonColors hcolors_;

    //store selected vertex
    int constraint_vhs_;


//    TriMesh& mesh_;
};
#endif
