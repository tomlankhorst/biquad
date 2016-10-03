#ifndef BIQUAD_BIQUAD_H
#define BIQUAD_BIQUAD_H

#include <vector>
#include <complex>

class BiQuadChain;

/** BiQuad class implements a single filter
 *
 * author: T.J.W. Lankhorst <t.j.w.lankhorst@student.utwente.nl>
 *
 * Filters that - in the z domain - are the ratio of two quadratic functions. The general form is:
 *
 *        b0 + b1 z^-1 + b2 z^-2
 * H(z) = ----------------------
 *        a0 + a1 z^-1 + a2 z^-2
 *
 * Which is often normalized by dividing all coefficients by a0.
 *
 * Example:
 * @code
 * #include "mbed.h"
 * #include <complex>
 *
 * // Example: 4th order Butterworth LP (w_c = 0.1*f_nyquist)
 * BiQuad bq1( 4.16599e-04, 8.33198e-04, 4.16599e-04, -1.47967e+00, 5.55822e-01 );
 * BiQuad bq2( 1.00000e+00, 2.00000e+00, 1.00000e+00, -1.70096e+00, 7.88500e-01 );
 *
 * BiQuadChain bqc;
 *
 * int main() {
 *
 *    // Add the biquads to the chain
 *    bqc.add( &bq1 ).add( &bq2 );
 *
 *    // Find the poles of the filter
 *    std::cout << "Filter poles" << std::endl;
 *    std::vector< std::complex<double> > poles = bqc.poles();
 *    for( size_t i = 0; i < poles.size(); i++ )
 *        std::cout << "\t"  << poles[i] << std::endl;
 *
 *    // Find the zeros of the filter
 *    std::cout << "Filter zeros" << std::endl;
 *    std::vector< std::complex<double> > zeros = bqc.zeros();
 *    for( size_t i = 0; i < poles.size(); i++ )
 *        std::cout << "\t" << zeros[i] << std::endl;
 *
 *    // Is the filter stable?
 *    std::cout << "This filter is " << (bqc.stable() ? "stable" : "instable") << std::endl;
 *
 *    // Output the step-response of 20 samples
 *  std::cout << "Step response 20 samples" << std::endl;
 *  for( int i = 0; i < 20; i++ )
 *      std::cout << "\t" << bqc.step( 1.0 ) << std::endl;
 * }
 * @endcode
 *
 * https://github.com/tomlankhorst/biquad
 *
 */
class BiQuad {

private:

    double B[3];
    double A[2];
    double wz[2];

    bool resetStateOnGainChange;

    /**
     * Sets the gain parameters
     */
    void set( double b0, double b1, double b2, double a1, double a2 );

public:

    /**
     * Initialize a unity TF biquad
     * @return BiQuad instance
     */
    BiQuad( );

    /**
     * Initialize a normalized biquad filter
     * @param b0
     * @param b1
     * @param b2
     * @param a1
     * @param a2
     * @return BiQuad instance
     */
    BiQuad( double b0, double b1, double b2, double a1, double a2 );

    /**
     * Initialize a biquad filter with all six coefficients
     * @param b0
     * @param b1
     * @param b2
     * @param a0
     * @param a1
     * @param a2
     * @return BiQuad instance
     */
    BiQuad( double b0, double b1, double b2, double a0, double a1, double a2 );

    /**
     * Initialize a PIDF biquad.
     * Based on Tustin-approx (trapezoidal) of the continous time version.
     * Behaviour equivalent to the PID controller created with the following MATLAB expression:
     *
     * C = pid( Kp, Ki, Kd, 1/N, Ts, 'IFormula', 'Trapezoidal', 'DFormula', 'Trapezoidal' );
     *
     * @param Kp    Proportional gain
     * @param Ki    Integral gain
     * @param Kd    Derivative gain
     * @param N     Filter coefficient ( N = 1/Tf )
     * @param Ts    Timestep
     */
    void PIDF( double Kp, double Ki, double Kd, double N, double Ts  );

    /**
     * Execute one digital timestep and return the result...
     * @param x input of the filer
     * @return output of the filter
     */
    double step( double x );

    /**
     * Return poles of the BiQuad filter
     * @return vector of std::complex poles
     */
    std::vector< std::complex<double> > poles( );

    /**
     * Return zeros of the BiQuad filter
     * @return vector of std::complex zeros
     */
    std::vector< std::complex<double> > zeros( );

    /**
     * Is this biquad stable?
     * Checks if all poles lie within the unit-circle
     * @return boolean whether the filter is stable or not
     */
    bool stable ();

    /**
     * Determines if the state variables are reset to zero on gain change.
     * Can be used for changing gain parameters on the fly.
     * @param v Value of the reset boolean
     */
    void setResetStateOnGainChange( bool v );

};

/**
 * The BiQuadChain class implements a chain of BiQuad filters
 */
class BiQuadChain {

private:
    std::vector< BiQuad* > biquads;
    std::vector< std::complex<double> > poles_zeros( bool zeros = false );

public:

    /**
     * Add a BiQuad pointer to the list: bqc.add(&bq);
     * @param bq Pointer to BiQuad instance
     * @return Pointer to BiQuadChain
     */
    BiQuadChain &add( BiQuad *bq );

    /**
     * Execute a digital time step cascaded through all bq's
     * @param x Input of the filter chain
     * @return Output of the chain
     */
    double step(double x);

    /**
     * Return poles of the BiQuad filter
     * @return vector of std::complex poles
     */
    std::vector< std::complex<double> > poles( );

    /**
     * Return zeros of the BiQuad filter
     * @return vector of std::complex zeros
     */
    std::vector< std::complex<double> > zeros( );

    /**
     * Is this biquad-chain stable?
     * Checks if all poles lie within the unit-circle
     * @return boolean whether the chain is stable or not
     */
    bool stable ();

    /**
     * Appends a BiQuad to the chain
     * Shorthand for .add(&bq)
     * @param bq BiQuad
     * @return Pointer to BiQuadChain
     */
    BiQuadChain &operator*( BiQuad& bq );

};

/**
 * Multiply two BiQuads
 * ... which in fact means appending them into a BiQuadChain
 * @return BiQuadChain of the two BiQuads
 */
BiQuadChain operator*( BiQuad&, BiQuad& );

#endif //BIQUAD_BIQUAD_H