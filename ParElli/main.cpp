#include "widget.h"
#include <QApplication>

int main(int argc, char *argv[])
{
    QApplication a(argc, argv);
    a.setWindowIcon(QIcon(QPixmap(":/icon/icon.png")));
    Widget w;
    w.show();

    w.updatePlotScale();

    return a.exec();
}
