#pragma once

#if QT_VERSION >= 0x050000
#include <QtWidgets>
#else
#include <QtGui>
#endif

#include "ui_MasterToolbarBase.h"


class MasterToolbar : public QWidget, public Ui::MasterToolbarBase
{
    Q_OBJECT

public:
    MasterToolbar(QWidget * parent = 0);
};