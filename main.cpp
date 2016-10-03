/**
 * Demo program for BiQuad and BiQuadChain classes
 * author: T.J.W. Lankhorst <t.j.w.lankhorst@student.utwente.nl>
 */

#include <stdlib.h>
#include <iostream>
#include <iomanip>
#include <complex>
#include "BiQuad.h"

// Example: 4th order Butterworth LP (w_c = 0.1*f_nyquist)
BiQuadChain bqc;
BiQuad bq1( 7.43110e-11, 1.64117e-10, 9.06793e-11, -1.45603e+00, 5.30959e-01 );
BiQuad bq2( 1.00000e+00, 2.14439e+00, 1.15577e+00, -1.47967e+00, 5.55816e-01 );
BiQuad bq3( 1.00000e+00, 2.04241e+00, 1.05319e+00, -1.52761e+00, 6.06223e-01 );
BiQuad bq4( 1.00000e+00, 1.93723e+00, 9.47394e-01, -1.60095e+00, 6.83335e-01 );
BiQuad bq5( 1.00000e+00, 1.85555e+00, 8.65227e-01, -1.70096e+00, 7.88499e-01 );
BiQuad bq6( 1.00000e+00, 1.81190e+00, 8.21307e-01, -1.82837e+00, 9.22458e-01 );

int main()
{

    bqc = bq1 * bq2 * bq3 * bq4 * bq5 * bq6;

    // Find the poles of the filter
    std::cout << "Filter poles" << std::endl;
    std::vector< std::complex<double> > poles = bqc.poles();
    for( size_t i = 0; i < poles.size(); i++ )
        std::cout << "\t"  << poles[i] << std::endl;

    // Find the zeros of the filter
    std::cout << "Filter zeros" << std::endl;
    std::vector< std::complex<double> > zeros = bqc.zeros();
    for( size_t i = 0; i < poles.size(); i++ )
        std::cout << "\t" << zeros[i] << std::endl;

    // Is the filter stable?
    std::cout << "This filter is " << (bqc.stable() ? "stable" : "instable") << std::endl;

    // Output the step-response of 20 samples
    std::cout << "Step response 20 samples" << std::endl;
    for( int i = 0; i < 20; i++ )
        std::cout << "\t" << bqc.step( 1.0 ) << std::endl;

    // Done
    return EXIT_SUCCESS;

}