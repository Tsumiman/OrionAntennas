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
    this->buildHyperbolaFocus();
    this->buildRayToParabolaFocus();
    this->buildHyperbola();
    this->buildNormalToHyperbola();
    this->buildHyperbolaCenter();
    this->buildRayToHyperbolaFocus();
    this->updateHyperbolaData();
    this->updateEfficiencyData();
    this->ui->plot->replot();
}

void Widget::updateHyperbolaData()
{
    double A = this->calculateHyperbolaParameterA();
    double C = this->calculateHyperbolaParameterC();
    double B = sqrt(C*C - A*A);
    double center = this->calculateHyperbolaCenter();
    this->ui->spinHA->setValue(A);
    this->ui->spinHB->setValue(B);
    this->ui->spinHC->setValue(C);
    this->ui->spinHE->setValue(C/A);
    this->ui->spinHCnt->setValue(center);

    double f = this->ui->spinf->value();
    QPointF intersectionPoint = findHyperbolaRayIntersection();
    double tgalpha = intersectionPoint.y() / (intersectionPoint.x() - f);
    this->ui->spinHAngle->setValue(atan(tgalpha) * 180 / M_PI * 2);
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
    double F = ui->spinF->value();
    int steps = 1000, current = 0;
    double step = D/steps;
    QVector<double> X(steps+1), Y(steps+1);
    for(double y = -D/2; y < D/2; y += step, ++current) {
        Y[current] = y;
        X[current] = y*y/4/F;
    }

    if (X.back() == 0 && Y.back() == 0) {
        X.pop_back();
        Y.pop_back();
    }
    this->Parabola->setData(Y, X);
}

void Widget::buildParabolaFocus()
{
    double F = ui->spinF->value();
    QVector<double> X(1), Y(1);
    X[0] = F; Y[0] = 0;
    this->Focus_Parabola->setData(Y, X);
}

void Widget::buildHyperbola()
{
    double a = calculateHyperbolaParameterA();
    double asq = a*a;
    double c = calculateHyperbolaParameterC();
    double csq = c*c;
    double bsq = csq - asq;
    double coef = asq/bsq;
    double center = (ui->spinF->value() + ui->spinf->value())/2;

    double d = ui->spind->value();
    int steps = 1000, current = 0;
    double step = d/steps;

    QVector<double> X(steps+1), Y(steps+1), shadowX(steps+1);
    for(double y = -d/2; y < d/2; y += step, ++current) {
        Y[current] = y;
        X[current] = sqrt(asq + coef*y*y) + center;
        shadowX[current] = -sqrt(asq+coef*y*y) + center;
    }

    if (X.back() == 0 && Y.back() == 0) {
        X.pop_back();
        Y.pop_back();
        shadowX.pop_back();
    }

    Hyperbola->setData(Y, X);
    Shadow_Hyperbola->setData(Y, shadowX);
}

void Widget::buildHyperbolaFocus()
{
    double f = ui->spinf->value();
    QVector<double> X(1), Y(1);
    X[0] = f; Y[0] = 0;
    this->Focus_Hyperbola->setData(Y, X);
}

void Widget::buildHyperbolaCenter()
{
    double center = this->calculateHyperbolaCenter();
    QVector<double> X(1), Y(1);
    X[0] = center; Y[0] = 0;
    this->Center_Hyperbola->setData(Y, X);
}

void Widget::buildRayToParabolaFocus()
{
    double D = ui->spinD->value();
    double F = ui->spinF->value();
    QVector<double> X(2), Y(2);
    QPointF intersection = this->findHyperbolaRayIntersection();
    X[0] = D*D/16/F, X[1] = intersection.x();
    Y[0] = D/2;      Y[1] = intersection.y();
    this->Ray_To_Parabola_Focus->setData(Y, X);
}

void Widget::buildRayToHyperbolaFocus()
{
    double f = ui->spinf->value();
    QPointF intersection = this->findHyperbolaRayIntersection();
    QVector<double> X(2), Y(2);
    X[0] = intersection.x(); X[1] = f;
    Y[0] = intersection.y(); Y[1] = 0;
    this->Reflected_Ray_To_Hyperbola_Focus->setData(Y, X);
}

void Widget::buildHyperbolaRayIntersection()
{
    QVector<double> X(1), Y(1);
    QPointF intersection = findHyperbolaRayIntersection();
    X[0] = intersection.x();
    Y[0] = intersection.y();
    this->Hyperbola_Ray_Intersection->setData(Y, X);
}

void Widget::buildNormalToHyperbola() {
    QPointF intersectionPoint = findHyperbolaRayIntersection();
    double a = calculateHyperbolaParameterA();
    double c = calculateHyperbolaParameterC();
    double center = (ui->spinF->value() + ui->spinf->value())/2;
    double x0 = intersectionPoint.x() - center;
    double y0 = intersectionPoint.y();

    double asq = a*a;
    double csq = c*c;
    double bsq = csq - asq;
    double coef = asq/bsq;

    double x1 = x0 - 1;
    double y1 = y0 - coef*y0/x0*(x1 - x0);

    QVector<double> X(2), Y(2);
    X[0] = x0+center, X[1] = x1+center;
    Y[0] = y0; Y[1] = y1;
    this->Normal_To_Hyperbola->setData(Y, X);
}

// ********************
// CALCULATIONS
// ********************

double Widget::calculateParabolaLength()
{
    double F = this->ui->spinF->value();
    double k = 1.0/4/F/F;
    double a = this->ui->spinD->value()/2;

    return (a + sqrt(1 + k*a*a) + log(sqrt(k)*a + sqrt(k*a*a + 1))/sqrt(k))/2;
}

double Widget::calculateIncommingArea()
{
    double D = this->ui->spinD->value();
    return (D/2) * (D/2) * M_PI;
}

QPointF Widget::findHyperbolaRayIntersection()
{
    double D = ui->spinD->value();
    double d = ui->spind->value();
    double F = ui->spinF->value();
    double tau = 1 - d/D;
    return QPointF((1-tau)*D*D/16/F + tau*F, d/2);
}

double Widget::calculateHyperbolaParameterA() {
    QPointF intersectionPoint = findHyperbolaRayIntersection();
    double center = (ui->spinF->value() + ui->spinf->value())/2;
    double c = calculateHyperbolaParameterC();
    double x0 = intersectionPoint.x() - center;
    double y0 = intersectionPoint.y();

    double A = 1;
    double B = (x0*x0 + y0*y0 + c*c);
    double C = x0*x0*c*c;
    double D = B*B - 4*A*C;

//    double a1 = sqrt((B + sqrt(D))/2);
    double a2 = sqrt((B - sqrt(D))/2);
    return a2;
}

double Widget::calculateHyperbolaParameterC()
{
    return (ui->spinF->value() - ui->spinf->value())/2;
}

double Widget::calculateHyperbolaCenter()
{
    return (ui->spinF->value() + ui->spinf->value())/2;
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

    this->initializeGraph(this->Hyperbola);
    this->Hyperbola->setName("Гиперболическое зеркало");
    this->initializeGraph(this->Center_Hyperbola);
    this->Center_Hyperbola->setName("Центр координат для гиперболы");
    this->initializeGraph(this->Shadow_Hyperbola);
    this->Shadow_Hyperbola->setName("Вторая ветка гиперболы");
    this->initializeGraph(this->Focus_Hyperbola);
    this->Focus_Hyperbola->setName("Фокус гиперболы");
    this->initializeGraph(this->Hyperbola_Ray_Intersection);
    this->Hyperbola_Ray_Intersection->setName("Точка пересечения луча с гиперболой");
    this->initializeGraph(this->Normal_To_Hyperbola);
    this->Normal_To_Hyperbola->setName("Нормаль в точке касания");
    this->initializeGraph(this->Reflected_Ray_To_Hyperbola_Focus);
    this->Reflected_Ray_To_Hyperbola_Focus->setName("Отраженный гиперболой луч");

    this->Parabola->setPen(QPen(QBrush(Qt::darkBlue), 2));
    this->Hyperbola->setPen(QPen(QBrush(Qt::darkRed), 2));
    this->Shadow_Hyperbola->setPen(QPen(QBrush(Qt::gray), 0.5));

    this->Focus_Parabola->setScatterStyle(QCPScatterStyle::ssCircle);
    this->Focus_Parabola->setPen(QPen(QBrush(Qt::red), 2));

    this->Hyperbola_Ray_Intersection->setScatterStyle(QCPScatterStyle::ssCircle);
    this->Hyperbola_Ray_Intersection->setPen(QPen(QBrush(Qt::black), 1));

    this->Center_Hyperbola->setScatterStyle(QCPScatterStyle::ssCircle);
    this->Center_Hyperbola->setPen(QPen(QBrush(Qt::black), 1));

    this->Focus_Hyperbola->setScatterStyle(QCPScatterStyle::ssCircle);
    this->Focus_Hyperbola->setPen(QPen(QBrush(Qt::darkGreen), 2));

    this->Ray_To_Parabola_Focus->setPen(QPen(QBrush(Qt::darkYellow), 0));
    this->Reflected_Ray_To_Hyperbola_Focus->setPen(QPen(QBrush(Qt::darkYellow), 0));
    this->Normal_To_Hyperbola->setPen(QPen(QBrush(Qt::darkCyan), 0));
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
    this->updatePlot();
}

void Widget::on_spinD_valueChanged(double)
{
    this->updatePlot();
}

void Widget::on_spinf_valueChanged(double)
{
    this->updatePlot();
}

void Widget::on_spind_valueChanged(double)
{
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
