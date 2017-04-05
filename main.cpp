/**
 * Demo program for BiQuad and BiQuadChain classes
 * author: T.J.W. Lankhorst <t.j.w.lankhorst@student.utwente.nl>
 */

#include <stdlib.h>
#include <iostream>
#include <iomanip>
#include <complex>
#include "BiQuad.h"

// Example: 3th order Butterworth LP (w_c = 0.1*f_nyquist)
BiQuadChain bqc;
BiQuad bq1( 3.40538e-04, 6.83088e-04, 3.42555e-04, -1.03207e+00, 2.75708e-01 );
BiQuad bq2( 1.00000e+00, 1.99997e+00, 9.99981e-01, -1.14298e+00, 4.12802e-01 );
BiQuad bq3( 1.00000e+00, 1.99412e+00, 9.94131e-01, -1.40438e+00, 7.35915e-01 );

int main()
{

    bqc = bq1 * bq2 * bq3;

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

    // Done0
    return EXIT_SUCCESS;

}