#ifndef _CONVECTION_SOLVER_
#define _CONVECTION_SOLVER_

#include <omp.h>
#include <complex>

#include "toefl/toefl.h"
#include "file/read_input.h"

enum target{
    TEMPERATURE, VORTICITY, POTENTIAL
};

struct Parameter
{
    double R; 
    double P; 
    double nu;
    double amp;
    unsigned nx; 
    unsigned nz; 
    double lx; 
    double lz; 
    double h; 
    double dt;
    unsigned itstp;
    enum toefl::bc bc_z;

    template< class Ostream>
    void display( Ostream& os)
    {
        os << "R is: "<< R <<"\n";
        os << "P is: "<< P <<"\n";
        os << "nu is: "<<nu  <<"\n";
        os << "nx is: "<< nx <<"\n";
        os << "nz is: "<< nz <<"\n";
        os << "lx is: "<< lx <<"\n";
        os << "lz is: "<< lz <<"\n";
        os << "h is: "<< h <<"\n";
        os << "dt is: "<< dt <<"\n";
    }
};

Parameter read( char const * file);

typedef std::complex<double> Complex;

void rayleigh_equations( toefl::QuadMat< Complex,2>& coeff, const Complex dx, const Complex dy, const Parameter& p);

inline void laplace_inverse( double& l_inv, const Complex dx, const Complex dy);
/*! @brief Solver for periodic boundary conditions of the toefl equations.
 * @ingroup solvers
 */
class Convection_Solver
{
  public:
    typedef toefl::Matrix<double, toefl::TL_DFT> Matrix_Type;
    /*! @brief Construct a solver for periodic boundary conditions
     *
     * The constructor allocates storage for the solver
     * and initializes all fourier coefficients as well as 
     * all low level solvers needed.  
     * @param blueprint Contains all the necessary parameters.
     * @throw Message If your parameters are inconsistent.
     */
    Convection_Solver( const Parameter& param);
    /*! @brief Prepare Solver for execution
     *
     * This function takes the fields and computes the missing 
     * one according to the target parameter passed. After that
     * it performs three initializing steps (one onestep-, 
     * one twostep-method and the threestep-method used in the step function)
     * in order to initialize the karniadakis scheme. The actual time is
     * thus T_0 + 3*dt after initialisation. 
     * @param v Container with three non void matrices
     * @param t which Matrix is missing?
     */
    void init( std::array< toefl::Matrix<double,toefl::TL_DFT>, 2>& v, enum target t);
    /*! @brief Perform a step by the 3 step Karniadakis scheme*/
    void step();
    /*! @brief Get the result
        
        You get the solution matrix of the current timestep.
        @param t The field you want
        @return A Read only reference to the field
        @attention The reference is only valid until the next call to 
            the step() function!
    */
    const Matrix_Type& getField( enum target t) const;
    /*! @brief Get the result

        Use this function when you want to call step() without 
        destroying the solution. 
        @param m 
            In exchange for the solution matrix you have to provide
            storage for further calculations. The field is swapped in.
        @param t 
            The field you want. 
        @attention The fields you get are not the ones of the current
            timestep. You get the fields that are not needed any more. 
            This means the densities are 4 timesteps "old" whereas 
            the potential is the one of the last timestep.
    */
    void getField( Matrix_Type& m, enum target t);
    /*! @brief Get the parameters of the solver.

        @return The parameters in use. 
        @note You cannot change parameters once constructed.
     */
    const Parameter& parameter() const { return param;}
    void setHeat( double x0, double y0, double sigma_x, double sigma_y, double amp);
  private:
    typedef std::complex<double> complex;
    //methods
    void init_coefficients( );
    void compute_cphi();//multiply cphi
    void first_steps(); 
    template< enum toefl::stepper S>
    void step_();
    //members
    const size_t rows, cols;
    const size_t crows, ccols;
    double x0_, y0_, sigma_x_, sigma_y_, amp_;
    const Parameter param;
    /////////////////fields//////////////////////////////////
    //GhostMatrix<double, TL_DFT> ghostdens, ghostphi;
    std::array< Matrix_Type, 2> dens, nonlinear;
    toefl::GhostMatrix<double, toefl::TL_DFT> phi;
    /////////////////Complex (void) Matrices for fourier transforms///////////
    std::array< toefl::Matrix< complex>, 2> cdens;
    toefl::Matrix< complex> cphi;
    ///////////////////Solvers////////////////////////
    toefl::Arakawa arakawa;
    toefl::Karniadakis<2, complex, toefl::TL_DFT> karniadakis;
    toefl::DFT_DRT dft_drt;
    /////////////////////Coefficients//////////////////////
    toefl::Matrix< double > phi_coeff;
};

#include "convection_solver.cpp"

#endif //_CONVECTION_SOLVER_
