#include <iostream>
#include "mainwindow.h"
#include <QApplication>

int main(int argc, char *argv[])
{
    QApplication a(argc, argv);
    MainWindow w; //is a derivative of QMainwindow & QWidget
    w.show(); //makes window visible
    //signal/slots get executed immediately, events in a chain
    return a.exec(); //enter qt's event loop
}
