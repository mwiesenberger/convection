#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QTimer>
#include <QMainWindow>
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
    cv::Mat matOriginal;
    cv::Mat matProcessed;
    cv::Mat matField;

    std::vector<cv::Vec3f> vecBlueCircles;
    std::vector<cv::Vec3f>::iterator itrBlueCircles;
    std::vector<cv::Vec3f> vecRedCircles;
    std::vector<cv::Vec3f>::iterator itrRedCircles;

    QTimer* tmrTimer;

public slots:
    void processFrameAndUpdateGUI();

private slots:
    void startButton_clicked();
    void resetButton_clicked();

};

#endif // MAINWINDOW_H
