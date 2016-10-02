#ifndef BIQUAD_BIQUAD_H
#define BIQUAD_BIQUAD_H

#include <vector>
#include <complex>

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
 * #include "BiQuad.h"
 *
 * BiQuadChain bqc;
 * BiQuad pidf;
 *
 * int main() {
 *     // Create a biquad filter based on PIDF parameters
 *    pidf.PIDF(1,1,1,1,1);
 *
 *    // Add the biquads to the chain
 *    bqc.add( &pidf );
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
     * Initialize a PIFD biquad.
     * Based on Tustin-approx (trapezoidal). of the continous time version
     * @param Kp
     * @param Ki
     * @param Kd
     * @param N
     * @param Ts
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
};

#endif //BIQUAD_BIQUAD_H