#include "window.h"
#include "ui_window.h"

Window::Window( QWidget *parent):
    QWidget( parent),
    ui( new Ui::Window)
{
    ui->setupUi(this);

}

Window::~Window()
{
    delete ui;
}

void Window::updateBuffer(cv::Mat next)
{
    QImage qImgCurrent( (uchar*)next.data, next.cols, next.rows, next.step, QImage::Format_RGB888);
    ui->label->setPixmap( QPixmap::fromImage( qImgCurrent));

}
