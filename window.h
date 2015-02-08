#ifndef WINDOW_H
#define WINDOW_H

#include <QWidget>
#include "opencv2/opencv.hpp"

namespace Ui{
    class Window;
}

class Window : public QWidget
{
    Q_OBJECT;
public:
    explicit Window(QWidget *parent);
    ~Window();
private:
    Ui::Window *ui;
public slots:
    void updateBuffer(cv::Mat next);
};

#endif // WINDOW_H
