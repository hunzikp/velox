#include <Rcpp.h>
#include <vector>
#include <math.h>

using namespace Rcpp;
using namespace std;

// [[Rcpp::export]]
NumericVector checktype_cpp(List rasterbands) {
  // out[0]: is integer
  // out[1]: has negative
  // out[2]: max abs value
  NumericVector out(3);
  out[0] = 1;
  double maxval = 0;
  int counter = 0;
  int nbands = rasterbands.size();
  for (int i = 0; i < nbands; i++) {
    NumericMatrix thisband = rasterbands[i];
    int size = thisband.nrow()*thisband.ncol();
    for (int j = 0; j < size; j++) {
      double val = thisband[j];
      if (fmod(val, 1) != 0) {
        out[0] = 0;
      }
      if ((counter == 0) || (abs(val) > maxval)) {
        maxval = abs(val);
      }
      if (val < 0) {
        out[1] = 1;
      }
      counter++;
    }
  }
  out[2] = maxval;
  return out;
}
