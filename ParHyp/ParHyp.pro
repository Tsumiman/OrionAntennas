#-------------------------------------------------
#
# Project created by QtCreator 2018-07-04T10:49:10
#
#-------------------------------------------------

QT       += core gui

greaterThan(QT_MAJOR_VERSION, 4): QT += widgets

TARGET = ParHyp
TEMPLATE = app
DEFINES += QT_NO_PRINTER

SOURCES += main.cpp\
        widget.cpp \
    qcustomplot.cpp

HEADERS  += widget.h \
    qcustomplot.h

FORMS    += widget.ui

RESOURCES += \
    icon.qrc
