#ifndef COLORS_H
#define COLORS_H

#include <QDialog>
#include <QTimer>
#include "opencv2/opencv.hpp"

namespace Ui {
class Colors;
}

class Colors : public QDialog
{
    Q_OBJECT

public:
    explicit Colors(QWidget *parent = 0);
    void getRanges( cv::Scalar& bmin, cv::Scalar& bmax);
    ~Colors();

private:
    Ui::Colors *ui;
    QTimer* timer;
    cv::VideoCapture webCam;
    cv::Mat current, currentHSV;
    cv::Mat red, blue;
    cv::Scalar blueMin, blueMax;

private slots:
    void processFrameAndUpdate();

};

#endif // COLORS_H
