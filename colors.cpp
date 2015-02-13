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
    blueMin = cv::Scalar( 110, 40 , 40);
    blueMax = cv::Scalar( 130, 255, 255);

}

Colors::~Colors()
{
    delete ui;
}

void Colors::getRanges(cv::Scalar &bmin, cv::Scalar &bmax)
{
    bmin = blueMin, bmax = blueMax;
}

void Colors::processFrameAndUpdate( )
{

    int w = ui->camera->width();
    int h = ui->camera->height();
    webCam >> current;
    cv::Mat shiftedHSV;
    cv::cvtColor(current, currentHSV, CV_BGR2HSV);
    cv::cvtColor(current, current, CV_BGR2RGB);
    cv::cvtColor(current, shiftedHSV, CV_BGR2HSV);
    //cv::cvtColor(current, current, CV_BGR2RGB);



    int bhmin = ui->bhmin->value(), bhmax = ui->bhmax->value();
    //setup dilate and erode and 60degree angle shift for red

    blueMin = cv::Scalar( bhmin, ui->bsmin->value(), ui->bvmin->value());
    blueMax = cv::Scalar( bhmax, ui->bsmax->value(), ui->bvmax->value());
    //redMin = cv::Scalar( ui->rhmin->value(), ui->rsmin->value(), ui->rvmin->value());
    //redMax = cv::Scalar( ui->rhmax->value(), ui->rsmax->value(), ui->rvmax->value());
    cv::inRange(currentHSV, blueMin, blueMax, blue );
    cv::inRange(shiftedHSV, blueMin, blueMax, red );

    int kkk = ui->dilute->value();
    cv::Size ksize( 2*kkk+1, 2*kkk+1 );
    cv::Point anchor( kkk, kkk );
    cv::Mat element = cv::getStructuringElement(cv::MORPH_ELLIPSE, ksize, anchor);

    cv::erode(  red, red, element, anchor, 1 );
    //ksize = cv::Size( 2*kkk + 9, 2*kkk+9);
    //ksanchor.x = ( kkk+4);
    //ksanchor.y = ( kkk+4);
    //element = cv::getStructuringElement(cv::MORPH_ELLIPSE, ksize, anchor);
    cv::dilate( red, red, element, anchor, 1 );

    cv::erode(  blue, blue, element, anchor, 1);
    cv::dilate( blue, blue, element, anchor, 1);
    QPixmap qPixCurrent = QPixmap::fromImage(QImage((uchar*)current.data, current.cols, current.rows, current.step, QImage::Format_RGB888 ));
    ui->camera->setPixmap(qPixCurrent.scaled(w,h,Qt::KeepAspectRatio));
    QPixmap qPixRed = QPixmap::fromImage(QImage((uchar*)red.data, red.cols, red.rows, red.step, QImage::Format_Indexed8));
    ui->redImage->setPixmap(qPixRed.scaled(w, h, Qt::KeepAspectRatio));
    QPixmap qPixBlue = QPixmap::fromImage(QImage((uchar*)blue.data, blue.cols, blue.rows, blue.step, QImage::Format_Indexed8));
    ui->blueImage->setPixmap(qPixBlue.scaled(w,h,Qt::KeepAspectRatio));

}
