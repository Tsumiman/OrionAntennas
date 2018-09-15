#include "widget.h"
#include "ui_widget.h"
#include <QDebug>
#include <QRect>
#include <QWheelEvent>

Widget::Widget(QWidget *parent) :
    QWidget(parent),
    ui(new Ui::Widget)
{
    ui->setupUi(this);

    this->initializeGraphs();
    this->ui->plot->setAttribute(Qt::WA_TransparentForMouseEvents);
    updatePlotScale();
    updatePlot();

//    this->ui->spinf->setMaximum(this->ui->spinF->value());
//    this->ui->spinF->setMinimum(this->ui->spinf->value());
    this->ui->spinD->setMinimum(this->ui->spind->value());
    this->ui->spind->setMaximum(this->ui->spinD->value());
}

Widget::~Widget()
{
    delete ui;
}

// ********************
// OUTPUT UPDATING
// ********************

void Widget::updatePlot()
{
    this->buildParabola();
    this->buildParabolaFocus();
    this->buildEllipsSecondFocus();
    this->buildRayToParabolaFocus();
    this->buildEllips();
    this->buildEllipsCenter();
    this->updateEllipsData();
    this->updateEfficiencyData();
    this->ui->plot->replot();
}

double Widget::calculateEllipsRotationAngle()
{
    double F = this->ui->spinF->value();
    double f = this->ui->spinf->value();
    double d = this->ui->spind->value();
    double tanalpha = (F-f)/(d/2);
    return atan(tanalpha);
}

double Widget::calculateEllipsParameterA()
{
    QPointF F1 = this->parabolaFocusPoint();
    QPointF F2 = this->ellipsSecondFocusPoint();
    QPointF M  = this->findEllipsRayIntersection();
    double a = (QLineF(F1, M).length() + QLineF(F2, M).length())/2;
    return a;
}

double Widget::calculateEllipsParameterB()
{
    double a = this->calculateEllipsParameterA();
    double c = this->calculateEllipsParameterC();
    return sqrt(a*a - c*c);
}

double Widget::calculateEllipsParameterC()
{
    QPointF F1 = this->parabolaFocusPoint();
    QPointF F2 = this->ellipsSecondFocusPoint();
    double c = (QLineF(F1, F2).length())/2;
    return c;
}

void Widget::updateEllipsData()
{
    QPointF center = this->calculateEllipsCenter();
    this->ui->spinECntX->setValue(center.x());
    this->ui->spinECntY->setValue(center.y());

    double a = this->calculateEllipsParameterA();
    this->ui->spinEA->setValue(a);

    double b = this->calculateEllipsParameterB();
    this->ui->spinEB->setValue(b);

    double c = this->calculateEllipsParameterC();
    this->ui->spinEC->setValue(c);

    double e = c/a;
    this->ui->spinEE->setValue(e);

    double alpha = this->calculateEllipsRotationAngle();
    this->ui->spinEAxisAngle->setValue(alpha/M_PI*180);
    double A = cos(alpha);
    double B = sin(alpha);

    // Ellips equation: x^2/b^2 + y^2/b^2 = 1
    // But x and y are rotated by angle alpha
    // and translated by vector "center".
    // Thus, to obtain x from y, we have to
    // thanslate and rotate back.
    double cA = (A*A/b/b + B*B/a/a);
    double cB = 2*A*B*((1/a/a) - (1/b/b));
    double cC = (A*A/a/a + B*B/b/b);
    double d = this->ui->spind->value();
    double u = d/2 - center.y();
    double D = u*u*cB*cB - 4*cA*(u*u*cC - 1);
    double z = (sqrt(D) - u*cB)/2/cA;
    double x = z + center.x();
    double f = this->ui->spinf->value();
    double angle = atan(d/2/(x-f));
    this->ui->spinEAngle->setValue(angle/M_PI*180);
}

void Widget::updateEfficiencyData()
{
    double frequency = this->ui->spinR->value();
    double c = 299.792458;
    double L = c/frequency;
    this->ui->spinL->setValue(L);

    double S = this->calculateIncommingArea();
    this->ui->spinS->setValue(S);

    double d = this->ui->spind->value();
    double D = this->ui->spinD->value();
    double A1 = 1.0;
    double alpha = d/D;
    double xi = (2*alpha*alpha/(2-A1))*(1 - alpha*alpha*A1/2);
    double zyu = (1 - xi); zyu = zyu*zyu;
    double q = this->ui->spinq->value();
    double G = 10*log10(4*M_PI/L/L*S*zyu*q);
    this->ui->spinG->setValue(G);

    this->ui->spinDL->setValue(d/L);
}

// ********************
// BUILDING GRAPHS
// ********************

void Widget::buildParabola()
{
    double D = ui->spinD->value();
    double d = ui->spind->value();
    double F = ui->spinF->value();
    int steps = 500, current = 0;
    double step = (D-d)/2/steps;
    QVector<double> X(steps+1), Y(steps+1);
    for(double y = d/2; current < steps; y += step, ++current) {
        Y[current] = -y;
        X[current] = (y-d/2)*(y-d/2)/4/F;
    }

    for (int i = 0; i < steps; ++i) {
        X.push_back(X[i]);
        Y.push_back(-Y[i]);
    }

    this->Parabola->setData(Y, X);
}

void Widget::buildParabolaFocus()
{
    QPointF F = this->parabolaFocusPoint();
    QVector<double> X(1), Y(1);
    X[0] = F.x(); Y[0] = F.y();
    this->Focus_Parabola->setData(Y, X);
}

void Widget::buildEllips()
{
    double a = this->calculateEllipsParameterA();
    double b = this->calculateEllipsParameterB();
    double alpha = this->calculateEllipsRotationAngle();
    double A = cos(alpha);
    double B = sin(alpha);
    QPointF center = this->calculateEllipsCenter();

    double cA = (A*A/b/b + B*B/a/a);
    double cB = 2*A*B*((1/a/a) - (1/b/b));
    double cC = (A*A/a/a + B*B/b/b);

    int steps = 500, current = 0;
    double d = ui->spind->value();
    double step = d/2/steps;
    QVector<double> X(500), Y(500);

    // Ellips equation: x^2/b^2 + y^2/b^2 = 1
    // But x and y are rotated by angle alpha
    // and translated by vector "center".
    // Thus, to obtain x from y, we have to
    // thanslate and rotate back.
    for (double y = 0; current < steps; y += step, ++current) {
        double u = y - center.y();
        double D = u*u*cB*cB - 4*cA*(u*u*cC - 1);
        double z = (sqrt(D) - u*cB)/2/cA;
        X[current] = z + center.x();
        Y[current] = -y;
    }
    for (int i = 0; i < steps; ++i) {
        X.push_back(X[i]);
        Y.push_back(-Y[i]);
    }

    this->Ellips->setData(Y, X);
}

void Widget::buildEllipsSecondFocus()
{
    QPointF F = this->ellipsSecondFocusPoint();
    QVector<double> X(1), Y(1);
    X[0] = F.x(); Y[0] = F.y();
    this->Focus_Ellips->setData(Y, X);
}

void Widget::buildEllipsCenter()
{
    QPointF center = this->calculateEllipsCenter();
    QVector<double> X(1), Y(1);
    X[0] = center.x(); Y[0] = center.y();
    this->Center_Ellips->setData(Y, X);
}

void Widget::buildRayToParabolaFocus()
{
    QPointF F1 = this->parabolaEndpoint();
    QPointF M = this->findEllipsRayIntersection();
    QVector<double> X(2), Y(2);
    X[0] = F1.x(), X[1] = M.x();
    Y[0] = F1.y(); Y[1] = M.y();
    this->Ray_To_Parabola_Focus->setData(Y, X);
}

// ********************
// CALCULATIONS
// ********************

double Widget::calculateIncommingArea()
{
    double D = this->ui->spinD->value();
    return (D/2) * (D/2) * M_PI;
}

QPointF Widget::findEllipsRayIntersection()
{
    QPointF F1 = this->parabolaFocusPoint();
    QPointF F2 = this->ellipsSecondFocusPoint();
    QPointF E = this->parabolaEndpoint();
    QPointF Z(-10,0);
    QPointF result;
    QLineF(E, F1).intersect(QLineF(Z, F2), &result);
    return result;
}

QPointF Widget::parabolaFocusPoint()
{
    double d = ui->spind->value();
    double F = ui->spinF->value();
    return QPointF(F, d/2);
}

QPointF Widget::ellipsSecondFocusPoint()
{
    double f = ui->spinf->value();
    return QPointF(f, 0);
}

QPointF Widget::calculateEllipsCenter()
{
    QPointF F1 = this->parabolaFocusPoint();
    QPointF F2 = this->ellipsSecondFocusPoint();
    return (F1+F2)/2;
}

QPointF Widget::parabolaEndpoint()
{
    double D = ui->spinD->value();
    double d = ui->spind->value();
    double F = ui->spinF->value();
    return QPointF((D-d)*(D-d)/16/F, D/2);
}

// ********************
// PLOTTING WIDGET SETUP
// ********************

void Widget::initializeGraphs() {
    this->initializeGraph(this->Parabola);
    this->Parabola->setName("Параболическое зеркало");
    this->initializeGraph(this->Focus_Parabola);
    this->Focus_Parabola->setName("Фокус параболы");
    this->initializeGraph(this->Ray_To_Parabola_Focus);
    this->Ray_To_Parabola_Focus->setName("Луч к фокусу параболы");

    this->initializeGraph(this->Ellips);
    this->Ellips->setName("Гиперболическое зеркало");
    this->initializeGraph(this->Center_Ellips);
    this->Center_Ellips->setName("Центр координат для гиперболы");
    this->initializeGraph(this->Focus_Ellips);
    this->Focus_Ellips->setName("Фокус гиперболы");

    this->Parabola->setPen(QPen(QBrush(Qt::darkBlue), 2));
    this->Ellips->setPen(QPen(QBrush(Qt::darkRed), 2));

    this->Focus_Parabola->setScatterStyle(QCPScatterStyle::ssCircle);
    this->Focus_Parabola->setPen(QPen(QBrush(Qt::red), 2));

    this->Center_Ellips->setScatterStyle(QCPScatterStyle::ssCircle);
    this->Center_Ellips->setPen(QPen(QBrush(Qt::black), 1));

    this->Focus_Ellips->setScatterStyle(QCPScatterStyle::ssCircle);
    this->Focus_Ellips->setPen(QPen(QBrush(Qt::darkGreen), 2));

    this->Ray_To_Parabola_Focus->setPen(QPen(QBrush(Qt::darkYellow), 0));
}

void Widget::initializeGraph(QCPGraph*& graph) {
    auto xAxis = this->ui->plot->xAxis;
    auto yAxis = this->ui->plot->yAxis;
    graph = new QCPGraph(xAxis, yAxis);
}

void Widget::updatePlotScale()
{
    double w = this->ui->plot->width();
    double h = this->ui->plot->height();
    this->ui->plot->xAxis->setRange(-w/this->xPerWidth/2, w/this->xPerWidth/2);
    this->ui->plot->yAxis->setRange(-h/this->yPerHeight/6, 5*h/this->yPerHeight/6);
    this->ui->plot->replot();
}

// ********************
// EVENTS
// ********************

void Widget::resizeEvent(QResizeEvent *event)
{
    QWidget::resizeEvent(event);
    this->updatePlotScale();
}

void Widget::wheelEvent(QWheelEvent *event)
{
    QWidget::wheelEvent(event);
    QPoint mousePosition = event->pos();
    QPoint widgetPosition = findPosition(this->ui->plot);
    QSize widgetSize = this->ui->plot->size();
    QRect rect(widgetPosition, widgetSize);
    if (!rect.contains(mousePosition)) {
        return;
    }
    if (event->angleDelta().y() > 0) {
        this->xPerWidth *= 1.1;
        this->yPerHeight *= 1.1;
    }
    else {
        this->xPerWidth /= 1.1;
        this->yPerHeight /= 1.1;
    }
    this->updatePlotScale();
}

QPoint Widget::findPosition(QWidget *child)
{
    if (child->parentWidget()) {
        return child->pos() + findPosition(child->parentWidget());
    }
    else {
        return QPoint(0,0);
    }
}

// ********************
// SIGNALS
// ********************

void Widget::on_spinF_valueChanged(double)
{
    //this->ui->spinf->setMaximum(arg1);
    this->updatePlot();
}

void Widget::on_spinD_valueChanged(double arg1)
{
    this->ui->spind->setMaximum(arg1);
    this->updatePlot();
}

void Widget::on_spinf_valueChanged(double)
{
    //this->ui->spinF->setMinimum(arg1);
    this->updatePlot();
}

void Widget::on_spind_valueChanged(double arg1)
{
    this->ui->spinD->setMinimum(arg1);
    this->updatePlot();
}

void Widget::on_spinq_valueChanged(double)
{
    this->updateEfficiencyData();
}

void Widget::on_spinR_valueChanged(double)
{
    this->updateEfficiencyData();
}
