/**
 * Demo program for BiQuad and BiQuadChain classes
 * author: T.J.W. Lankhorst <t.j.w.lankhorst@student.utwente.nl>
 */

#include <stdlib.h>
#include <iostream>
#include <iomanip>
#include <complex>
#include "BiQuad.h"

int main()
{

    BiQuadChain bqc;

    // Filter consists of two biquad sections
    BiQuad bq1( 4.82434e-03, 9.64869e-03, 4.82434e-03, -1.04860e+00, 2.96140e-01 );
    BiQuad bq2( 1.00000e+00, 2.00000e+00, 1.00000e+00, -1.32091e+00, 6.32739e-01 );

    // Add the biquads to the chain
    bqc.add( &bq1 ).add( &bq2 );

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