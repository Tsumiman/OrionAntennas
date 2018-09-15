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


    QCPGraph *Parabola, *Hyperbola, *Focus_Parabola, *Center_Hyperbola,
     *Focus_Hyperbola, *Ray_To_Parabola_Focus, *Hyperbola_Ray_Intersection,
    *Normal_To_Hyperbola, *Reflected_Ray_To_Hyperbola_Focus, *Shadow_Hyperbola;

    double xPerWidth = 10;
    double yPerHeight = 10;

public:
    explicit Widget(QWidget *parent = 0);
    void updatePlotScale();
    ~Widget();

private:
    Ui::Widget *ui;

    void updatePlot();
    QPointF findHyperbolaRayIntersection();
    double calculateHyperbolaParameterA();
    double calculateHyperbolaParameterC();
    double calculateHyperbolaCenter();

    void buildParabola();
    void buildParabolaFocus();
    void buildHyperbola();
    void buildHyperbolaFocus();
    void buildHyperbolaCenter();
    void buildRayToParabolaFocus();
    void buildRayToHyperbolaFocus();
    void buildHyperbolaRayIntersection();
    void buildNormalToHyperbola();
    void updateHyperbolaData();
    void updateEfficiencyData();

    double calculateParabolaLength();
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
