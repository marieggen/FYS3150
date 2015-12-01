TEMPLATE = app
CONFIG += console c++11
CONFIG -= app_bundle
CONFIG -= qt

SOURCES += main.cpp

INCLUDEPATH += /usr/local/include # header files
LIBS += -L/usr/local/lib # library files
LIBS += -larmadillo -llapack -lblas
