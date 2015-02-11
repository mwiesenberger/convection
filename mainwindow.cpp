#include <QString>
//#include <omp.h>
#include "mainwindow.h"
#include "ui_mainwindow.h"
#include "convection_solver.h"

typedef typename Convection_Solver::Matrix_Type Matrix;

MainWindow::MainWindow(QWidget *parent) :
    QMainWindow(parent), window(new Window(0)),
    ui(new Ui::MainWindow), capWebCam(0)
{
    ui->setupUi(this);
    //capWebCam.open(0);
    if( !capWebCam.isOpened()){
        ui->textField->append("error: Webcam not accessed succesfully!\n");
        return;
    }
    capWebCam.set(CV_CAP_PROP_FRAME_WIDTH, 320);
    capWebCam.set(CV_CAP_PROP_FRAME_HEIGHT, 240);
    tmrTimer = new QTimer(this); //if this is deleted so will be Timer
    //default parameters
    Parameter p;
    p.P = 1;
    p.R = 500000;
    p.nu = 0.01;
    p.amp = 100000;
    p.nz = 32;
    p.nx = 256;
    p.dt = 4e-6;
    p.itstp = 5;
    p.bc_z = toefl::TL_DST10; //dst 2
    omp_set_num_threads( 1);
    p.lz = 1.;
    p.h = p.lz / (double)p.nz;
    p.lx = (double)p.nx * p.h;
    ui->textField->setText( "\
\n\
        * Input-File for CONVECTION - Code *\n\
        ------------------------------------\n\
\n\
1) pr (Prandtl number)                 = 1.0\n\
2) ra (rayleigh number)                = 500000\n\
3) nu (artif. hyper-viscosity)         = 0.01\n\
4) amp (initial perturbation amplitude)= 100000\n\
5) nz (grid points in z, best 2^n)     = 32\n\
6) nx (grid points in x, best 2^n)     = 256\n\
7) dt (time step)                      = 4e-6\n\
8) itstp (steps between output)        = 5\n\
9) OMP_NUM_THREADS                     = 1\n"
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
}

void MainWindow::processFrameAndUpdateGUI(){
    capWebCam >> matOriginal;
    if( matOriginal.empty()) return;

    cv::inRange( matOriginal, cv::Scalar( 0,0,175), cv::Scalar( 100, 100, 256), matProcessed); //if Pixel is in Range 1, 0 else, type CV_8U, (tracks red color)
    cv::GaussianBlur( matProcessed, matProcessed, cv::Size(9,9), 1.5); //smooth edges
    cv::HoughCircles( matProcessed, vecRedCircles, CV_HOUGH_GRADIENT, 2, matProcessed.rows/4, 200, 100, 0,0); //vecRedCircles contains x,y,r for each detected circle, minimum distance between circles = rows/4, min/max_radius to be detected = 0 (unknown)
    for( itrRedCircles=vecRedCircles.begin(); itrRedCircles!=vecRedCircles.end(); itrRedCircles++)
    {
        ui->textField->append(
                  QString("ball position: x = ") + (QString::number( (*itrRedCircles)[0])).rightJustified(4, ' ')
                + QString(" y = ") + (QString::number( (*itrRedCircles)[1])).rightJustified(4, ' ')
                + QString(" ball radius: = ") + (QString::number( (*itrRedCircles)[2], 'f', 3)).rightJustified(7, ' '));

        cv::circle(matOriginal, cv::Point((int)(*itrRedCircles)[0], (int)(*itrRedCircles)[1]), 3, cv::Scalar( 0, 255, 0), -1);
        cv::circle(matOriginal, cv::Point((int)(*itrRedCircles)[0], (int)(*itrRedCircles)[1]), (int)(*itrRedCircles)[2] , cv::Scalar( 0, 0, 255), 1);
    }
    cv::cvtColor( matOriginal, matOriginal, CV_BGR2RGB); //change for qt rgb image data type

    QImage qImgOriginal( (uchar*)matOriginal.data, matOriginal.cols, matOriginal.rows, matOriginal.step, QImage::Format_RGB888);
    ui->webcam->setPixmap( QPixmap::fromImage(qImgOriginal));


    const Matrix* field = &solver->getField( TEMPERATURE);
    std::vector<double> visual( field->rows()*field->cols());

    //add gradient and swap z direction
    for( unsigned i=0; i<field->rows(); i++)
        for( unsigned j=0; j<field->cols(); j++)
            visual[i*field->cols()+j] = (*field)(field->rows() -  1 - i,j) + solver->parameter().R*(-0.5+(double)i/(double)(field->rows() + 1));

    window->updateBuffer( visual, field->cols(), field->rows(), solver->parameter().R/2.);

    for(unsigned i=0; i<solver->parameter().itstp; i++)
    {
        solver->step(); //here is the timestep
    }


}

void MainWindow::startButton_clicked() {
    if( !tmrTimer->isActive())
    {
        tmrTimer->start(20); //starts the timer and emit a timeout signal every 20ms (50Hz), but only when event loop (exec()) is running
        window->show();
        ui->startButton->setText("pause");
    }
    else
    {
        tmrTimer->stop();
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
