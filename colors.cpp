#include "colors.h"
#include "ui_colors.h"

Colors::Colors(QWidget *parent) :
    QDialog(parent),
    ui(new Ui::Colors)
{
    webCam.open(0);
    webCam.set(CV_CAP_PROP_FRAME_WIDTH, 320);
    webCam.set(CV_CAP_PROP_FRAME_HEIGHT, 240);
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

void Colors::getRanges(cv::Scalar &bmin, cv::Scalar &bmax, int& dilute)
{
    bmin = blueMin, bmax = blueMax, dilute = kkk;
}

void Colors::processFrameAndUpdate( )
{

    int w = ui->camera->width();
    int h = ui->camera->height();
    webCam >> current;
    cv::Mat shiftedHSV;
    cv::cvtColor(current, currentHSV, CV_BGR2HSV);
    cv::cvtColor(current, current,    CV_BGR2RGB);
    cv::cvtColor(current, shiftedHSV, CV_BGR2HSV);

    int bhmin = ui->bhmin->value(), bhmax = ui->bhmax->value();

    blueMin = cv::Scalar( bhmin, ui->bsmin->value(), ui->bvmin->value());
    blueMax = cv::Scalar( bhmax, ui->bsmax->value(), ui->bvmax->value());
    //redMin = cv::Scalar( ui->rhmin->value(), ui->rsmin->value(), ui->rvmin->value());
    //redMax = cv::Scalar( ui->rhmax->value(), ui->rsmax->value(), ui->rvmax->value());
    cv::inRange(currentHSV, blueMin, blueMax, blue );
    cv::inRange(shiftedHSV, blueMin, blueMax, red );
    //erode and dilate to reduce noise
    kkk = ui->dilute->value();
    cv::Size ksize( 2*kkk+1, 2*kkk+1 );
    cv::Point anchor( kkk, kkk );
    cv::Mat element = cv::getStructuringElement(cv::MORPH_ELLIPSE, ksize, anchor);
    cv::erode(  red, red, element, anchor, 1 );
    cv::dilate( red, red, element, anchor, 1 );
    cv::erode(  blue, blue, element, anchor, 1);
    cv::dilate( blue, blue, element, anchor, 1);
    //find circles
    std::vector<cv::Vec3f> redC, blueC;
    cv::GaussianBlur( red, red, cv::Size(7,7), 0,0); //smooth edges
    cv::GaussianBlur( blue, blue, cv::Size(7,7), 0,0); //smooth edges
    cv::HoughCircles( red, redC, CV_HOUGH_GRADIENT, 2, red.rows/4, 200, 100, 0,0); //vecRedCircles contains x,y,r for each detected circle, minimum distance between circles = rows/4, min/max_radius to be detected = 0 (unknown)
    cv::HoughCircles( blue, blueC, CV_HOUGH_GRADIENT, 2, blue.rows/4, 200, 100, 0,0); //vecRedCircles contains x,y,r for each detected circle, minimum distance between circles = rows/4, min/max_radius to be detected = 0 (unknown)
    //show circles
    for( unsigned i=0; i<redC.size(); i++)
    {
        cv::circle(current, cv::Point((int)redC[i][0], (int)redC[i][1]), 7, cv::Scalar( 0, 0, 255), -1);
        cv::circle(current, cv::Point((int)redC[i][0], (int)redC[i][1]), (int)redC[i][2] , cv::Scalar( 0, 0, 255), 1);
    }

    for( unsigned i=0; i<blueC.size(); i++)
    {
        cv::circle(current, cv::Point((int)blueC[i][0], (int)blueC[i][1]), 7, cv::Scalar( 255, 0, 0), -1);
        cv::circle(current, cv::Point((int)blueC[i][0], (int)blueC[i][1]), (int)blueC[i][2] , cv::Scalar( 255, 0, 0), 1);
    }
    //show images
    QPixmap qPixCurrent = QPixmap::fromImage(QImage((uchar*)current.data, current.cols, current.rows, current.step, QImage::Format_RGB888 ));
    ui->camera->setPixmap(qPixCurrent.scaled(w,h,Qt::KeepAspectRatio));
    QPixmap qPixRed = QPixmap::fromImage(QImage((uchar*)red.data, red.cols, red.rows, red.step, QImage::Format_Indexed8));
    ui->redImage->setPixmap(qPixRed.scaled(w, h, Qt::KeepAspectRatio));
    QPixmap qPixBlue = QPixmap::fromImage(QImage((uchar*)blue.data, blue.cols, blue.rows, blue.step, QImage::Format_Indexed8));
    ui->blueImage->setPixmap(qPixBlue.scaled(w,h,Qt::KeepAspectRatio));

}
