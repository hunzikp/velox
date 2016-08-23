#include <Rcpp.h>
#include <vector>
#include "median.h"

using namespace Rcpp;
using namespace std;


// [[Rcpp::export]]
NumericMatrix medianfocal_cpp(NumericMatrix rasterband, int wrow, int wcol, int band) {

  int nrow = rasterband.nrow();
  int ncol = rasterband.ncol();

  NumericMatrix newband(nrow, ncol);
  vector<double> v;

  int idim = (wrow-1)/2;
  int jdim = (wcol-1)/2;

  for (int i = 0; i < nrow; i++) {
    int imin = i - idim;
    int imax = i + idim + 1;
    for (int j = 0; j < ncol; j++) {
      int jmin = j - jdim;
      int jmax = j + jdim + 1;
      for (int wi = imin; wi < imax; wi++) {
        for (int wj = jmin; wj < jmax; wj++) {
          if (wi >= 0 && wi < nrow && wj >= 0 && wj < ncol) {
            v.push_back(rasterband(wi, wj));
          }
        }
      }
      newband(i,j) = median(v);
      v.clear();
    }
  }
  return(newband);
}

// [[Rcpp::export]]
NumericMatrix sumfocal_cpp(NumericMatrix rasterband, NumericMatrix weights, int wrow, int wcol, int band) {

  int ncol = rasterband.ncol();
  int nrow = rasterband.nrow();
  NumericMatrix newband(nrow, ncol);

  int idim = (wrow-1)/2;
  int jdim = (wcol-1)/2;

  for (int i = 0; i < nrow; i++) {
    int imin = i - idim;
    int imax = i + idim + 1;
    for (int j = 0; j < ncol; j++) {
      int jmin = j - jdim;
      int jmax = j + jdim + 1;
      double cusum = 0;
      int wicounter = 0;
      for (int wi = imin; wi < imax; wi++) {
        int wjcounter = 0;
        for (int wj = jmin; wj < jmax; wj++) {
          if (wi >= 0 && wi < nrow && wj >= 0 && wj < ncol) {
            cusum = cusum + rasterband(wi, wj)*weights(wicounter, wjcounter);
          }
          wjcounter++;
        }
        wicounter++;
      }
      newband(i,j) = cusum;
    }
  }
  return(newband);
}

// [[Rcpp::export]]
NumericMatrix meanfocal_cpp(NumericMatrix rasterband, NumericMatrix weights, int wrow, int wcol, int band) {

  int ncol = rasterband.ncol();
  int nrow = rasterband.nrow();
  NumericMatrix newband(nrow, ncol);

  int idim = (wrow-1)/2;
  int jdim = (wcol-1)/2;

  for (int i = 0; i < nrow; i++) {
    int imin = i - idim;
    int imax = i + idim + 1;
    for (int j = 0; j < ncol; j++) {
      int jmin = j - jdim;
      int jmax = j + jdim + 1;
      double cusum = 0;
      int cucount = 0;
      int wicounter = 0;
      for (int wi = imin; wi < imax; wi++) {
        int wjcounter = 0;
        for (int wj = jmin; wj < jmax; wj++) {
          if (wi >= 0 && wi < nrow && wj >= 0 && wj < ncol) {
            cusum = cusum + rasterband(wi, wj)*weights(wicounter, wjcounter);
            cucount++;
          }
          wjcounter++;
        }
        wicounter++;
      }
      newband(i,j) = cusum/cucount;
    }
  }
  return(newband);
}


