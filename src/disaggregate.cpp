#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix disaggregate_cpp(NumericMatrix band, IntegerVector factor) {
  NumericVector dim = band.attr("dim"); 
  
  int nrow = dim[0];
  int ncol = dim[1];
  
  int drow = factor[0];
  int dcol = factor[1];
  
  NumericMatrix newband = no_init(nrow*drow, ncol*dcol);
  
  for (int j = 0; j < ncol; j++) {
    for (int i = 0; i < nrow; i++) {
      for (int q = 0; q < dcol; q++) {
        for (int p = 0; p < drow; p++) {
          newband(i*drow + p, j*dcol + q) = band(i, j);
        }
      }        
    }
  }
  
  return newband;
}