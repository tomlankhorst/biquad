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
     */
    BiQuad( );
    /**
     * Initialize a normalized biquad filter
     */
    BiQuad( double b0, double b1, double b2, double a1, double a2 );
    /**
     * Initialize a biquad filter with all six coefficients
     */
    BiQuad( double b0, double b1, double b2, double a0, double a1, double a2 );
    /**
     * Initialize a PIFD biquad based on Tustin-approx (trapezoidal). of the continous time version
     */
    void PIDF( double Kp, double Ki, double Kd, double N, double Ts  );
    /**
     * Execute one digital timestep and return the result...
     */
    double step( double x );
    /**
     * Return poles of the BiQuad filter
     */
    std::vector< std::complex<double> > poles( );
    /**
     * Return zeros of the BiQuad filter
     */
    std::vector< std::complex<double> > zeros( );
    /**
     * Is this biquad stable? Checks if all poles lie within the unit-circle
     */
    bool stable ();
    /**
     * Determines if the state variables are reset to zero on gain change.
     * Can be used for changing gain parameters on the fly.
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
    /** Add a BiQuad pointer to the list: bqc.add(&bq);
     */
    BiQuadChain &add( BiQuad *bq );
    /** Execute a digital time step cascaded through all bq's
     */
    double step(double x);
    /** Return poles of the BiQuad filter
     */
    std::vector< std::complex<double> > poles( );
    /** Return zeros of the BiQuad filter
     */
    std::vector< std::complex<double> > zeros( );
    /** Is this biquad stable? Checks if all poles lie within the unit-circle
     */
    bool stable ();
};


#endif //BIQUAD_BIQUAD_H
