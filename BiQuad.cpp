#include "BiQuad.h"

BiQuad::BiQuad(double b0, double b1, double b2, double a1, double a2) {

    B[0] = b0; B[1] = b1; B[2] = b2;
    A[0] = a1; A[1] = a2;

    wz[0] = 0; wz[1] = 0;

}

BiQuad::BiQuad(double b0, double b1, double b2, double a0, double a1, double a2) {
    // Do the normalization...
    BiQuad(b0/a0,b1/a0,b2/a0,a1/a0,a2/a0);
}

double BiQuad::step(double x) {

    double y,w;

    /* Direct form II */
    w =      x - A[0]*wz[0] - A[1]*wz[1];
    y = B[0]*w + B[1]*wz[0] + B[2]*wz[1];

    /* Shift */
    wz[1] = wz[0];
    wz[0] = w;

    return y;

}

std::vector< std::complex<double> > BiQuad::poles() {

    std::vector< std::complex<double> > poles;

    std::complex<double> b2(A[0]*A[0],0);
    std::complex<double> ds = std::sqrt( b2-4*A[1] );

    poles.push_back( 0.5*(-A[0]+ds) );
    poles.push_back( 0.5*(-A[0]-ds) );

    return poles;

}

std::vector< std::complex<double> > BiQuad::zeros() {

    std::vector< std::complex<double> > zeros;

    std::complex<double> b2(B[1]*B[1],0);
    std::complex<double> ds = std::sqrt( b2-4*B[0]*B[2] );

    zeros.push_back( 0.5*(-B[1]+ds)/B[0] );
    zeros.push_back( 0.5*(-B[1]-ds)/B[0] );

    return zeros;

}

BiQuadChain &BiQuadChain::add(BiQuad *bq) {
    biquads.push_back( bq );
    return *this;
}

double BiQuadChain::step(double x) {

    int i;
    size_t bqs;

    bqs = biquads.size();

    for( i = 0; i < bqs; i++ )
        x = biquads[i]->step( x );

    return x;
}

std::vector< std::complex<double> > BiQuadChain::poles_zeros( bool zeros ) {

    std::vector< std::complex<double> > chain, bq;
    int i;
    size_t bqs;

    bqs = biquads.size();

    for( i = 0; i < bqs; i++ ){
        bq = ( zeros ) ? biquads[ i ]->zeros() : biquads[ i ]->poles();
        chain.insert( chain.end(), bq.begin(), bq.end() );
    }

    return chain;

}

std::vector< std::complex<double> > BiQuadChain::poles() {
    return poles_zeros( false );
}

std::vector< std::complex<double> > BiQuadChain::zeros() {
    return poles_zeros( true );
}