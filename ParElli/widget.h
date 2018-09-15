#ifndef WIDGET_H
#define WIDGET_H

#include <QWidget>
#include "qcustomplot.h"

namespace Ui {
class Widget;
}

class Widget : public QWidget
{
    Q_OBJECT


    QCPGraph *Parabola, *Ellips, *Focus_Parabola, *Center_Ellips,
     *Focus_Ellips, *Ray_To_Parabola_Focus;

    double xPerWidth = 10;
    double yPerHeight = 10;

public:
    explicit Widget(QWidget *parent = 0);
    void updatePlotScale();
    ~Widget();

private:
    Ui::Widget *ui;

    void updatePlot();

    double calculateEllipsRotationAngle();
    double calculateEllipsParameterA();
    double calculateEllipsParameterB();
    double calculateEllipsParameterC();

    QPointF findEllipsRayIntersection();
    QPointF parabolaFocusPoint();
    QPointF ellipsSecondFocusPoint();
    QPointF calculateEllipsCenter();
    QPointF parabolaEndpoint();

    double calculateHyperbolaParameterA();
    double calculateHyperbolaParameterC();

    void buildParabola();
    void buildParabolaFocus();
    void buildEllips();
    void buildEllipsSecondFocus();
    void buildEllipsCenter();
    void buildRayToParabolaFocus();
    void updateEllipsData();
    void updateEfficiencyData();

    double calculateIncommingArea();

    void initializeGraphs();
    void initializeGraph(QCPGraph *&graph);

protected:
    void resizeEvent(QResizeEvent *event) override;
    void wheelEvent(QWheelEvent *event) override;
    QPoint findPosition(QWidget *child);

private slots:
    void on_spinF_valueChanged(double arg1);
    void on_spinD_valueChanged(double arg1);
    void on_spinf_valueChanged(double arg1);
    void on_spind_valueChanged(double arg1);
    void on_spinq_valueChanged(double arg1);
    void on_spinR_valueChanged(double arg1);
};

#endif // WIDGET_H
