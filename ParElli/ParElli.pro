#-------------------------------------------------
#
# Project created by QtCreator 2018-07-04T10:49:10
#
#-------------------------------------------------

QT       += core gui printsupport

greaterThan(QT_MAJOR_VERSION, 4): QT += widgets

TARGET = ParElli
TEMPLATE = app

SOURCES += main.cpp\
        widget.cpp \
    qcustomplot.cpp

HEADERS  += widget.h \
    qcustomplot.h

FORMS    += widget.ui

RESOURCES += \
    icon.qrc

CONFIG += c++11
