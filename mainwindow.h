#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QTimer>
#include <QMainWindow>
#include <QFileDialog>
#include <QTextStream>
#include "window.h"
#include "opencv2/opencv.hpp"

//#include "convection_solver.h"

class Convection_Solver;
typedef Convection_Solver Solver;
//typedef typename Solver::Matrix_Type Matrix;

namespace Ui {
    class MainWindow; //forward declare the ui functionality class
}

class MainWindow : public QMainWindow
{
    Q_OBJECT

public:
    explicit MainWindow(QWidget *parent = 0);
    ~MainWindow();

private:
    Ui::MainWindow *ui; //gives access to ui widgets from designer via pimpl
    Window * window;

    Solver * solver;

    cv::VideoCapture capWebCam;
    cv::Mat current;
    cv::Mat matProcessed;
    cv::Mat matField;
    cv::Scalar blueMin, blueMax;

    cv::Mat red;
    cv::Mat blue;
    int dilute;

    QTimer* tmrTimer;

public slots:
    void processFrameAndUpdateGUI();

private slots:
    void startButton_clicked();
    void resetButton_clicked();
    void saveAction_triggered();
    void loadAction_triggered();
    void exitAction_triggered();

    void on_actionSetColors_triggered();
};

#endif // MAINWINDOW_H
