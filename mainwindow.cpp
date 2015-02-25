#include <QString>
#include <omp.h>
#include "mainwindow.h"
#include "ui_mainwindow.h"
#include "convection_solver.h"
#include "colors.h"

typedef typename Convection_Solver::Matrix_Type Matrix;

MainWindow::MainWindow(QWidget *parent) :
    QMainWindow(parent), window(new Window(0)),
    ui(new Ui::MainWindow), capWebCam(0)
{
    ui->setupUi(this); //creates user interface
    //capWebCam.open(0);
    if( !capWebCam.isOpened()){
        ui->textField->append("error: Webcam not accessed succesfully!\n");
        return;
    }
    capWebCam.set(CV_CAP_PROP_FRAME_WIDTH, 320);
    capWebCam.set(CV_CAP_PROP_FRAME_HEIGHT, 240);
    tmrTimer = new QTimer(this); //if this is deleted so will be Timer
    tmrTimer->start(20);
    blueMin = cv::Scalar( 110, 40 , 40);
    blueMax = cv::Scalar( 130, 255, 255);
    //default parameters
    Parameter p;
    p.P = 1;
    p.R = 5e5;
    p.nu = 1e-5;
    p.amp = 1e5;
    p.nz = 32;
    p.nx = 256;
    p.dt = 4e-6;
    p.itstp = 20;
    p.bc_z = toefl::TL_DST10; //dst 2
    omp_set_num_threads( 2);
    p.lz = 1.;
    p.h = p.lz / (double)p.nz;
    p.lx = (double)p.nx * p.h;
    ui->textField->setText( "\
\n\
        * Input-File for CONVECTION - Code *\n\
        ------------------------------------\n\
\n\
1) pr (Prandtl number)                 = 1.0\n\
2) ra (rayleigh number)                = 5e5\n\
3) nu (artif. hyper-viscosity)         = 1e-5\n\
4) amp (initial perturbation amplitude)= 1e5\n\
5) nz (grid points in z, best 2^n)     = 32\n\
6) nx (grid points in x, best 2^n)     = 256\n\
7) dt (time step)                      = 4e-6\n\
8) itstp (steps between output)        = 20\n\
9) OMP_NUM_THREADS                     = 2\n"
);

    solver = new Solver( p);

    //init to zero
    Matrix theta( p.nz, p.nx, 0.), vorticity( theta);
    // initialize theta here ...
    init_gaussian( theta, 0.2,0.5, 5./128./p.lx, 5./128., p.amp);
    init_gaussian( theta, 0.7,0.3, 5./128./p.lx, 5./128., p.amp);
    init_gaussian( theta, 0.9,0.2, 5./128./p.lx, 5./128., -p.amp);
    //initialize solver
    try{
        std::array< Matrix,2> arr{{ theta, vorticity}};
        //now set the field to be computed
        solver->init( arr, POTENTIAL);
    }catch( toefl::Message& m){m.display();}

    connect( tmrTimer, SIGNAL(timeout()), this, SLOT( processFrameAndUpdateGUI()));
    connect( ui->startButton, SIGNAL(clicked()), this, SLOT( startButton_clicked()));
    connect( ui->resetButton, SIGNAL(clicked()), this, SLOT( resetButton_clicked()));

    connect( ui->actionSave, SIGNAL(triggered()), this, SLOT(saveAction_triggered()));
    connect( ui->actionLoad, SIGNAL(triggered()), this, SLOT(loadAction_triggered()));
    connect( ui->actionExit, SIGNAL(triggered()), this, SLOT(exitAction_triggered()));

}

MainWindow::~MainWindow()
{
    delete ui;
    delete solver;
    delete window;
    delete tmrTimer;
}

void MainWindow::processFrameAndUpdateGUI(){
    capWebCam >> current;
    cv::Mat shiftedHSV, currentHSV;
    cv::cvtColor(current, currentHSV, CV_BGR2HSV);
    cv::cvtColor(current, current,    CV_BGR2RGB);
    cv::cvtColor(current, shiftedHSV, CV_BGR2HSV);
    cv::inRange(currentHSV, blueMin, blueMax, blue );
    cv::inRange(shiftedHSV, blueMin, blueMax, red );
    //erode and dilate to reduce noise
    int kkk = dilute;
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

    //show circles in scaled image
    unsigned nx = solver->parameter().nx;
    unsigned nz = solver->parameter().nz;
    unsigned nxC = current.cols;
    unsigned nzC = current.rows;
    double lx = solver->parameter().lx;
    double lz = solver->parameter().lz;

    cv::Mat share(nz, nx, CV_8UC3, cv::Scalar(0, 0, 0));
    //show circles
    for( unsigned i=0; i<redC.size(); i++)
        cv::circle(share, cv::Point((int)(redC[i][0]/nxC*nx), (int)(redC[i][1]/nzC*nz)), 4, cv::Scalar( 255, 0, 0), -1);

    for( unsigned i=0; i<blueC.size(); i++)
        cv::circle(share, cv::Point((int)(blueC[i][0]/nxC*nx), (int)(blueC[i][1]/nzC*nz)), 4, cv::Scalar( 0, 0, 255), -1);
    QPixmap pixmap = QPixmap::fromImage(QImage( share.data, nx, nz, QImage::Format_RGB888));
    int w = ui->webcam_2->width();
    int h = ui->webcam_2->height();
    ui->webcam_2->setPixmap( pixmap.scaled( w, h, Qt::KeepAspectRatio, Qt::SmoothTransformation));
    //////////////////////show webcam image and found circles
    int webCamWidth = current.cols;
    int webCamHeight = current.rows;
    QPixmap qPixCurrent = QPixmap::fromImage(QImage((uchar*)current.data, current.cols, current.rows, current.step, QImage::Format_RGB888 ));
    ui->webcam->setPixmap(qPixCurrent.scaled(webCamWidth,webCamHeight,Qt::KeepAspectRatio));

    //show simulation and make some more steps//////////

    if( window->isVisible())
    {
        const Matrix* field = &solver->getField( TEMPERATURE);
        std::vector<double> visual( field->rows()*field->cols());
        //add gradient and swap z direction
        for( unsigned i=0; i<field->rows(); i++)
            for( unsigned j=0; j<field->cols(); j++)
                visual[i*field->cols()+j] = (*field)(field->rows() -  1 - i,j) + solver->parameter().R*(-0.5+(double)i/(double)(field->rows() + 1));
        //show simulation
        window->updateBuffer( visual, field->cols(), field->rows(), solver->parameter().R/2.);

        Matrix src( solver->parameter().nz, solver->parameter().nx, 0.);
        for( unsigned i=0; i<redC.size(); i++)
            init_gaussian( src, 1.-redC[i][0]/nxC , 1.-redC[i][1]/nzC, 0.07/lx, 0.07, solver->parameter().amp);

        for( unsigned i=0; i<blueC.size(); i++)
            init_gaussian( src, 1.-blueC[i][0]/nxC, 1.-blueC[i][1]/nzC, 0.07/lx, 0.07, -solver->parameter().amp);

        for(unsigned i=0; i<solver->parameter().itstp; i++)
        {
            solver->step(src); //here is the timestep
        }
    }


}

void MainWindow::startButton_clicked() {
    if( !window->isVisible())
    {
        //tmrTimer->start(20); //starts the timer and emit a timeout signal every 20ms (50Hz), but only when event loop (exec()) is running
        window->show();
        ui->startButton->setText("pause");
    }
    else
    {
        window->hide();
        //tmrTimer->stop();
        ui->startButton->setText("resume");
    }
}

void MainWindow::resetButton_clicked() {
    window->hide();
    ui->startButton->setText("start");
    QString input = ui->textField->toPlainText();
    Parameter p = read( input.toStdString());
    //reallocate ressources
    delete solver;
    solver = new Solver( p);

    //init to zero
    Matrix theta( p.nz, p.nx, 0.), vorticity( theta);
    //initialize solver
    try{
        // initialize theta here ...
        init_gaussian( theta, 0.2,0.5, 5./128./p.lx, 5./128., p.amp);
        init_gaussian( theta, 0.7,0.3, 5./128./p.lx, 5./128., p.amp);
        init_gaussian( theta, 0.9,0.2, 5./128./p.lx, 5./128., -p.amp);
        std::array< Matrix,2> arr{{ theta, vorticity}};
        //now set the field to be computed
        solver->init( arr, POTENTIAL);
    }catch( toefl::Message& m){m.display();}
}

void MainWindow::saveAction_triggered()
{
    QString saveFile = QFileDialog::getSaveFileName(this);
    QString input = ui->textField->toPlainText();
    QFile data( saveFile);
    if( data.open( QFile::WriteOnly | QFile::Truncate))
    {
        QTextStream out(&data);
        out << input;
    }
}

void MainWindow::loadAction_triggered()
{
    QString loadFile = QFileDialog::getOpenFileName(this);
    std::string input = file::read_file( loadFile.toStdString().data());
    QString qInput(QString::fromStdString(input));
    ui->textField->setText(qInput);
}
void MainWindow::exitAction_triggered()
{
    window->close();
    this->close();
}

void MainWindow::on_actionSetColors_triggered()
{
    tmrTimer->stop();
    capWebCam.release();
    Colors* dialog = new Colors(this);
    dialog->setModal(true);
    dialog->exec();
    if( dialog->result() == QDialog::Accepted )
        dialog->getRanges(blueMin, blueMax, dilute);
    delete dialog;
    capWebCam.open(0);
    capWebCam.set(CV_CAP_PROP_FRAME_WIDTH, 320);
    capWebCam.set(CV_CAP_PROP_FRAME_HEIGHT, 240);
    tmrTimer->start(20);

}
