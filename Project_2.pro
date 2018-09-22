TEMPLATE = app
CONFIG += console c++11
CONFIG -= app_bundle
CONFIG -= qt

INCLUDEPATH += C:\Armadillo\include

SOURCES += \
        jacobi_matrix.cpp \

x64: LIBS += \
    -L\C:\Users\Elisabeth\Documents\cpp1\armadillo\examples\lib_win64\
    -lblas_win64_MT \
    -llapack_win64_MT

HEADERS += \
    catch2.h \
    offdiag.h \
    jacobi_method.h
