#include <Rcpp.h>

using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix getcoordinates_cpp(NumericVector dim, NumericVector res, NumericVector extent){

  int nrow = dim[0];
  int ncol = dim[1];
  double xres = res[0];
  double yres = res[1];
  double xmin = extent[0];
  double ymax = extent[3];

  NumericMatrix coordmat = NumericMatrix(nrow*ncol, 2);
  int counter = 0;
  for (int i=0; i<nrow; i++) {
    for (int j=0; j<ncol; j++) {
      coordmat(counter, 0) =  xmin + j*xres + (xres/2);
      coordmat(counter, 1) =  ymax - i*yres - (yres/2);
      counter++;
    }
  }
  return coordmat;
}
