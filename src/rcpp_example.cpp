// load Rcpp
#include <Rcpp.h>

using namespace Rcpp;		// shorthand

// [[Rcpp::export]]
NumericMatrix RcppGibbs(int n, int thn) {

    int i,j;
    NumericMatrix mat(n, 2);

    // The rest of the code follows the R version
    double x=0, y=0;

    for (i=0; i<n; i++) {
        for (j=0; j<thn; j++) {
            x = R::rgamma(3.0,1.0/(y*y+4));
            y = R::rnorm(1.0/(x+1),1.0/sqrt(2*x+2));
        }
        mat(i,0) = x;
        mat(i,1) = y;
    }

    return mat;             // Return to R
}
