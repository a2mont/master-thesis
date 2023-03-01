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
class MasterPlugin : public QObject, BaseInterface, ToolboxInterface, LoggingInterface
{
Q_OBJECT
    Q_INTERFACES(BaseInterface)
    Q_INTERFACES(ToolboxInterface)
    Q_INTERFACES(LoggingInterface)

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


private slots:
    // initialization functions
    void initializePlugin();

    void simpleLaplace();


public :
  
  ~MasterPlugin() {};
  
  QString name() { return QString("Master Plugin"); };
  
  QString description() { return QString("Plugin used for my thesis !"); };
  

private:
    QSpinBox* iterationsSpinbox_;
};
#endif