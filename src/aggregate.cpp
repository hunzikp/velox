#include <Rcpp.h>
#include <vector>
#include "median.h"

using namespace Rcpp;
using namespace std;

// [[Rcpp::export]]
List aggregate_cpp(List rasterbands, NumericVector dim, NumericVector res, NumericVector factor, int aggtype) {
  // aggtype 0: sum <<- default if aggtype outside allowed range
  // aggtype 1: mean
  // aggtype 2: min
  // aggtype 3: max
  // aggtype 4: median

  int nrow = dim[0];
  int ncol = dim[1];
  int nbands = rasterbands.size();

  int nrown = (int)floor(nrow/factor[1]);
  int ncoln = (int)floor(ncol/factor[0]);
  int startrow, endrow, startcol, endcol;
  Rcpp::List newrasterbands(nbands);

  // Outer raster band loop
  for (int k = 0; k < nbands; k++) {

    NumericMatrix thisband = rasterbands[k];
    NumericMatrix newband(nrown, ncoln);

    if (aggtype <= 1 || aggtype > 4) { // Sum-type aggregation
      double cusum, counter;
      for (int i = 0; i < nrown; i++) {
        startrow = (int)rint(i*factor[1]);
        endrow = (int)(rint((i+1)*factor[1]) - 1);
        for (int j = 0; j < ncoln; j++) {
          startcol = (int)rint(j*factor[0]);
          endcol = (int)(rint((j+1)*factor[0]) - 1);
          cusum = 0;
          counter = 0;
          for (int in = startrow; in <= endrow; in++) {
            for (int jn = startcol; jn <= endcol; jn++) {
              cusum = cusum + thisband(in, jn);
              counter++;
            }
          }
          if (aggtype == 1) {
            newband(i, j) = cusum/counter;
          } else {
            newband(i, j) = cusum;
          }
        }
      }
    } else if (aggtype == 2 || aggtype == 3) { // Extremal aggregation
      double val, minval, maxval;
      for (int i = 0; i < nrown; i++) {
        startrow = (int)rint(i*factor[1]);
        endrow = (int)(rint((i+1)*factor[1]) - 1);
        for (int j = 0; j < ncoln; j++) {
          startcol = (int)rint(j*factor[0]);
          endcol = (int)(rint((j+1)*factor[0]) - 1);
          minval = thisband(startrow, startcol);
          maxval = minval;
          for (int in = startrow; in <= endrow; in++) {
            for (int jn = startcol; jn <= endcol; jn++) {
              val = thisband(in, jn);
              if (val < minval) {
                minval = val;
              }
              if (val > maxval) {
                maxval = val;
              }
            }
          }
          if (aggtype == 2) {
            newband(i, j) = minval;
          } else {
            newband(i, j) = maxval;
          }
        }
      }
    } else if (aggtype == 4) { // Median aggregation
      vector<double> valvec;
      for (int i = 0; i < nrown; i++) {
        startrow = (int)rint(i*factor[1]);
        endrow = (int)(rint((i+1)*factor[1]) - 1);
        for (int j = 0; j < ncoln; j++) {
          startcol = (int)rint(j*factor[0]);
          endcol = (int)(rint((j+1)*factor[0]) - 1);
          for (int in = startrow; in <= endrow; in++) {
            for (int jn = startcol; jn <= endcol; jn++) {
              valvec.push_back(thisband(in, jn));
            }
          }
          newband(i, j) = median(valvec);
          valvec.clear();
        }
      }
    }
    // Add aggregated band to new raster band list
    newrasterbands[k] = newband;
  }

  double xres = res[0]*factor[0];
  double yres = res[1]*factor[1];
  NumericVector newres = NumericVector::create(xres, yres);
  NumericVector newdim = NumericVector::create(nrown, ncoln);

  Rcpp::List out(3);
  out[0] = newrasterbands;
  out[1] = newdim;
  out[2] = newres;

  return(out);
}
