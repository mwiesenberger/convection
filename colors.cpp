#include "colors.h"
#include "ui_colors.h"

Colors::Colors(QWidget *parent) :
    QDialog(parent),
    ui(new Ui::Colors)
{
    webCam.open(0);
    ui->setupUi(this);
    timer = new QTimer(this);
    timer->start(20);
    connect(timer, SIGNAL(timeout()), this, SLOT(processFrameAndUpdate()));
}

Colors::~Colors()
{
    delete ui;
}
void Colors::getRanges(cv::Scalar &rmin, cv::Scalar &rmax, cv::Scalar &bmin, cv::Scalar &bmax)
{
    rmin = redMin, rmax = redMax, bmin = blueMin, bmax = blueMax;
}

void Colors::processFrameAndUpdate( )
{
    int w = ui->camera->width();
    int h = ui->camera->height();
     webCam >> current;
     cv::cvtColor(current, currentHSV, CV_BGR2HSV);
     cv::cvtColor(current, current, CV_BGR2RGB);


     int bhmin = ui->bhmin->value(), bhmax = ui->bhmax->value();
     //setup dilate and erode and 60degree angle shift for red

     blueMin = cv::Scalar( bhmin, ui->bsmin->value(), ui->bvmin->value());
     blueMax = cv::Scalar( bhmax, ui->bsmax->value(), ui->bvmax->value());
     redMin = cv::Scalar( ui->rhmin->value(), ui->rsmin->value(), ui->rvmin->value());
     redMax = cv::Scalar( ui->rhmax->value(), ui->rsmax->value(), ui->rvmax->value());
     cv::inRange(currentHSV, blueMin, blueMax, blue );
     cv::inRange(currentHSV, redMin, redMax, red );

     QPixmap qPixCurrent = QPixmap::fromImage(QImage((uchar*)current.data, current.cols, current.rows, current.step, QImage::Format_RGB888 ));
     ui->camera->setPixmap(qPixCurrent.scaled(w,h,Qt::KeepAspectRatio));
     QPixmap qPixRed = QPixmap::fromImage(QImage((uchar*)red.data, red.cols, red.rows, red.step, QImage::Format_Indexed8));
     ui->redImage->setPixmap(qPixRed.scaled(w, h, Qt::KeepAspectRatio));
     QPixmap qPixBlue = QPixmap::fromImage(QImage((uchar*)blue.data, blue.cols, blue.rows, blue.step, QImage::Format_Indexed8));
     ui->blueImage->setPixmap(qPixBlue.scaled(w,h,Qt::KeepAspectRatio));


}
