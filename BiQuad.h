/**
 * BiQuad Filter Classes
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
 */

#ifndef BIQUAD_BIQUAD_H
#define BIQUAD_BIQUAD_H

#include <vector>
#include <complex>

/**
 * BiQuad class implements a single filter
 */
class BiQuad {

private:
    double B[3];
    double A[2];
    double wz[2];

public:
    /**
     * Initialize a normalized biquad filter
     */
    BiQuad( double b0, double b1, double b2, double a1, double a2 );
    /**
     * Initialize a biquad filter with all six coefficients
     */
    BiQuad( double b0, double b1, double b2, double a0, double a1, double a2 );
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
     */
    BiQuadChain &add( BiQuad *bq );
    /**
     * Execute a digital time step cascaded through all bq's
     */
    double step(double x);
    /**
     * Return poles of the BiQuad filter
     */
    std::vector< std::complex<double> > poles( );
    /**
     * Return zeros of the BiQuad filter
     */
    std::vector< std::complex<double> > zeros( );
};


#endif //BIQUAD_BIQUAD_H
