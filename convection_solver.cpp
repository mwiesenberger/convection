#include "convection_solver.h"
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif


Parameter read( const std::string& str)
{
    Parameter p;
    std::vector<double> para;
    try{ para = file::read_input( str); }
    catch (toefl::Message& m)
    {
        m.display();
        throw m;
    }
    p.P = para[1];
    p.R = para[2];
    p.nu = para[3];
    p.amp = para[4];
    p.nz = para[5];
    p.nx = para[6];
    p.dt = para[7];
    p.itstp = para[8];
    p.bc_z = toefl::TL_DST10;
    omp_set_num_threads( para[9]);
    //std::cout<< "With "<<omp_get_max_threads()<<" threads\n";

    p.lz = 1.;
    p.h = p.lz / (double)p.nz;
    p.lx = (double)p.nx * p.h;
    return p;
}
void rayleigh_equations( toefl::QuadMat< Complex,2>& coeff, const Complex dx, const Complex dy, const Parameter& p)
{
    double laplace = (dx*dx + dy*dy).real();
    coeff( 0,0) = laplace - p.nu*laplace*laplace,  coeff( 0,1) = -p.R*dx/laplace;
    coeff( 1,0) = -p.P*dx,    coeff( 1,1) = p.P*laplace - p.nu*laplace*laplace;
}
inline void laplace_inverse( double& l_inv, const Complex dx, const Complex dy)
{
    l_inv = 1.0/(dx*dx + dy*dy).real();
}

inline void Convection_Solver::step(){ step_<toefl::TL_ORDER3>();}
//void Convection_Solver::setHeat( std::vector<double> x0, std::vector<double> y0, std::vector<double> sigma_x, std::vector<double> sigma_y, std::vector<double> amp){
////if( x0 >= 0 && x0 < param.lx && y0 >= 0 && y0 < param.lz)
//    x0_ = x0, y0_ = y0, sigma_x_ = sigma_x, sigma_y_ = sigma_y, amp_ = amp;
////else
//    //std::cerr << "x0 or y0 is not within the boundaries!\n";
//}

Convection_Solver::Convection_Solver( const Parameter& p):
    rows( p.nz ), cols( p.nx ),
    crows( rows), ccols( cols/2+1),
    //x0_(0), amp_(0),
    param( p),
    //fields
    dens( toefl::MatrixArray<double, toefl::TL_DFT, 2>::construct( rows, cols)), nonlinear( dens),
    phi( rows, cols, param.bc_z, toefl::TL_PERIODIC),
    cdens( toefl::MatrixArray<complex, toefl::TL_NONE, 2>::construct( crows, ccols)),
    cphi(crows, ccols),
    //Solvers
    arakawa( p.h),
    karniadakis(rows, cols, crows, ccols, p.dt),
    dft_drt( rows, cols, fftw_convert( p.bc_z), FFTW_MEASURE),
    //Coefficients
    phi_coeff( crows, ccols)
{
    init_coefficients( );
}

void Convection_Solver::init_coefficients( )
{
    toefl::Matrix< toefl::QuadMat< complex, 2> > coeff( crows, ccols);
    const complex kxmin( 0, 2.*M_PI/param.lx), kzmin( 0, M_PI/param.lz);
    // dft_drt is not transposing so i is the y index by default
    for( unsigned i = 0; i<crows; i++)
        for( unsigned j = 0; j<ccols; j++)
        {
            rayleigh_equations( coeff( i,j), (double)j*kxmin, (double)(i+1)*kzmin, param);
            laplace_inverse( phi_coeff( i,j), (double)j*kxmin, (double)(i+1)*kzmin);
        }
    double norm = param.nx * fftw_normalisation( param.bc_z, param.nz);
    karniadakis.init_coeff( coeff, norm);
}
void Convection_Solver::init( std::array< toefl::Matrix<double, toefl::TL_DFT>,2>& v, enum target t)
{
    //fourier transform input into cdens
    for( unsigned k=0; k<2; k++)
    {
#ifdef TL_DEBUG
        if( v[k].isVoid())
            throw Message("You gave me a void Matrix!!", _ping_);
#endif
        dft_drt.r2c( v[k], cdens[k]);
    }
    //don't forget to normalize coefficients!!
    double norm = param.nx * fftw_normalisation( param.bc_z, param.nz);
    for( unsigned k=0; k<2; k++)
        for( unsigned i=0; i<crows; i++)
            for( unsigned j=0; j<ccols;j++)
                cdens[k](i,j) /= norm;
    switch( t) //which field must be computed?
    {
        case( TEMPERATURE):
            throw toefl::Message( "Temperature independent", _ping_);
            break;
        case( VORTICITY):
            //bring cdens and cphi in the right order
            swap_fields( cphi, cdens[1]);
            //solve for cdens[1]
            for( unsigned i=0; i<crows; i++)
                for( unsigned j=0; j<ccols; j++)
                    cdens[1](i,j) = cphi(i,j) /phi_coeff(i,j);
            break;
        case( POTENTIAL):
            //solve for cphi
            for( unsigned i=0; i<crows; i++)
                for( unsigned j=0; j<ccols; j++)
                {
                    cphi(i,j) = cdens[1](i,j)*phi_coeff(i,j);
                }
            break;
    }
    //backtransform to x-space
    for( unsigned k=0; k<2; k++)
    {
        dft_drt.c2r( cdens[k], dens[k]);
    }
    dft_drt.c2r( cphi, phi);
    //now the density and the potential is given in x-space
    first_steps();
}

void Convection_Solver::getField( toefl::Matrix<double, toefl::TL_DFT>& m, enum target t)
{
#ifdef TL_DEBUG
    if(m.isVoid())
        throw Message( "You may not swap in a void Matrix!\n", _ping_);
#endif
    switch( t)
    {
        case( TEMPERATURE):    swap_fields( m, nonlinear[0]); break;
        case( VORTICITY):      swap_fields( m, nonlinear[1]); break;
        case( POTENTIAL):      swap_fields( m, cphi); break;
    }
}
const toefl::Matrix<double, toefl::TL_DFT>& Convection_Solver::getField( enum target t) const
{
    toefl::Matrix<double, toefl::TL_DFT> const * m = 0;
    switch( t)
    {
        case( TEMPERATURE):     m = &dens[0]; break;
        case( VORTICITY):       m = &dens[1]; break;
        case( POTENTIAL):       m = &phi; break;
    }
    return *m;
}

void Convection_Solver::first_steps()
{
    karniadakis.invert_coeff<toefl::TL_EULER>( );
    step_<toefl::TL_EULER>();
    karniadakis.invert_coeff<toefl::TL_ORDER2>();
    step_<toefl::TL_ORDER2>();
    karniadakis.invert_coeff<toefl::TL_ORDER3>();
    step_<toefl::TL_ORDER3>();
}

void Convection_Solver::compute_cphi()
{
#pragma omp parallel for
    for( size_t i = 0; i < crows; i++)
        for( size_t j = 0; j < ccols; j++)
            cphi(i,j) = phi_coeff(i,j)*cdens[1](i,j);
}

inline void Convection_Solver::step(const Matrix_Type& src){ 
    step_<toefl::TL_ORDER3>(&src);
}
template< enum toefl::stepper S>
void Convection_Solver::step_(const Matrix_Type* src )
{
    phi.initGhostCells(  );
    //1. Compute nonlinearity
#pragma omp parallel for
    for( unsigned k=0; k<2; k++)
    {
        toefl::GhostMatrix<double, toefl::TL_DFT> ghostdens{ rows, cols, param.bc_z, toefl::TL_PERIODIC, toefl::TL_VOID};
        swap_fields( dens[k], ghostdens); //now dens[k] is void
        ghostdens.initGhostCells( );
        arakawa( phi, ghostdens, nonlinear[k]);
        swap_fields( dens[k], ghostdens); //now ghostdens is void
    }
    if( src != nullptr)
        for( unsigned i=0; i<rows; i++)
            for(unsigned j=0; j<cols; j++)
                dens[0](i,j) += src->operator()(i,j); //not exactly correct but who cares?
                //init_gaussian( dens[0], x0_[i][0], x0_[i][1], x0_[i][2], x0_[i][2], amp_[i]);
    //2. perform karniadakis step
    karniadakis.step_i<S>( dens, nonlinear);
    //3. solve linear equation
    //3.1. transform v_hut
#pragma omp parallel for
    for( unsigned k=0; k<2; k++)
        dft_drt.r2c( dens[k], cdens[k]);
    //3.2. perform karniadaksi step and multiply coefficients for phi
    karniadakis.step_ii( cdens);
    compute_cphi();
    //3.3. backtransform
#pragma omp parallel for
    for( unsigned k=0; k<2; k++)
        dft_drt.c2r( cdens[k], dens[k]);
    dft_drt.c2r( cphi,  phi); //field in phi again
}
