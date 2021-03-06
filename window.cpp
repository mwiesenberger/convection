
#include "window.h"
#include "ui_window.h"
#include "draw/colormap.h"

Window::Window( QWidget *parent):
    QWidget( parent),
    ui( new Ui::Window)
{
    ui->setupUi(this);
    map = new draw::ColorMapRedBlueExt;

}

Window::~Window()
{
    delete ui;
    delete map;
}

void Window::updateBuffer(std::vector<double>& show, unsigned width, unsigned height, double scale)
{
    map->scale() = scale;

    //convert according to color map
    std::vector<draw::Color> resource( show.size());
    for( unsigned i=0; i<show.size(); i++)
        resource[i] = (*map)( show[i]);
    cv::Mat share( height, width, CV_32FC3, (float*)resource.data());
    current.create( height, width, CV_8UC3);
    share.convertTo( current, CV_8UC3, 256.,0);

    QPixmap pixmap = QPixmap::fromImage(QImage( current.data, width, height, QImage::Format_RGB888));
    int w = ui->label->width();
    int h = ui->label->height();
    ui->label->setPixmap( pixmap.scaled( w, h, Qt::KeepAspectRatio, Qt::SmoothTransformation));
}
