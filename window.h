#ifndef WINDOW_H
#define WINDOW_H

#include <QWidget>
#include "opencv2/opencv.hpp"

namespace draw{
    class ColorMapRedBlueExt;}
namespace Ui{
    class Window;
}
class Convection_Solver;

class Window : public QWidget
{
    Q_OBJECT;
public:
    explicit Window(QWidget *parent);
    ~Window();
private:
    cv::Mat current;
    Ui::Window *ui;
    draw::ColorMapRedBlueExt * map;
    std::vector<double> visual;
public slots:
    void updateBuffer(std::vector<double> &show, unsigned width, unsigned height,double scale);
};

#endif // WINDOW_H
