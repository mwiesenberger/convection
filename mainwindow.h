#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QTimer>
#include <QMainWindow>
#include "opencv2/opencv.hpp"

//#include "convection_solver.h"

class Convection_Solver;
typedef Convection_Solver Solver;
//typedef typename Solver::Matrix_Type Matrix;

namespace Ui {
class MainWindow;
}

class MainWindow : public QMainWindow
{
    Q_OBJECT

public:
    explicit MainWindow(QWidget *parent = 0);
    ~MainWindow();

private:
    Ui::MainWindow *ui; //gives access to ui widgets from designer

    Solver * solver;

    cv::VideoCapture capWebCam;
    cv::Mat matOriginal;
    cv::Mat matProcessed;

    QImage qimageOriginal;
    QImage qimageProcessed;

    std::vector<cv::Vec3f> vecBlueCircles;
    std::vector<cv::Vec3f>::iterator itrBlueCircles;
    std::vector<cv::Vec3f> vecRedCircles;
    std::vector<cv::Vec3f>::iterator itrRedCircles;

    QTimer* tmrTimer;

public slots:
    void processFrameAndUpdateGUI();

private slots:
    void on_pushButton_clicked();
    void on_restart_clicked();

};

#endif // MAINWINDOW_H
