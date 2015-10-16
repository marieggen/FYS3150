TEMPLATE = app
CONFIG += console c++11
CONFIG -= app_bundle
CONFIG -= qt

SOURCES += main.cpp \
    functions.cpp \
    lib.cpp

INCLUDEPATH += /usr/local/include # header files
LIBS += -L/usr/local/lib # library files
# LIBS += -larmadillo -llapack -lblas

# Finner med kommandoen "mpic++ --showme:link"
INCLUDEPATH += -I/usr/local/Cellar/open-mpi/1.10.0/include
# Finner med kommandoen "mpic++ --showme:compile"
LIBS += -L/usr/local/opt/libevent/lib -L/usr/local/Cellar/open-mpi/1.10.0/lib -lmpi_cxx -lmpi

HEADERS += \
    functions.h \
    lib.h
