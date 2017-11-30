#include <Rcpp.h>
#include <vector>

using namespace Rcpp;
using namespace std;

// [[Rcpp::export]]
NumericMatrix pointextract_cpp(List rasterbands, NumericVector dim, NumericVector extent, NumericVector res, NumericMatrix pcoords) {

  double xmin = extent[0];
  double xmax = extent[1];
  double ymin = extent[2];
  double ymax = extent[3];

  int nbands = rasterbands.size();

  double xres = res[0];
  double yres = res[1];

  NumericMatrix out(pcoords.nrow(), nbands);

  for (int i = 0; i < pcoords.nrow(); i++) {
    double x = pcoords(i,0);
    double y = pcoords(i,1);
    if (x > xmin && xmax > x && y > ymin && ymax > y) {
      double rdist = ymax - y;
      int row = floor(rdist/yres);
      double cdist = x - xmin;
      int col = floor(cdist/xres);
      for (int k = 0; k < nbands; k++) {
        NumericMatrix this_band;
        this_band = as<NumericMatrix>(rasterbands[k]);
        out(i,k) = this_band(row, col);
      }
    } else {
      for (int k = 0; k < nbands; k++) {
        out(i,k) = NA_REAL;
      }
    }
  }

  return out;
}
