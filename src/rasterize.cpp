#include <Rcpp.h>

using namespace Rcpp;

NumericMatrix coord2index(NumericMatrix coordmat, NumericVector extent, NumericVector res) {

  double xmin = extent[0];
  double ymax = extent[3];
  double xres = res[0];
  double yres = res[1];

  NumericMatrix indexmat = NumericMatrix(coordmat.nrow(), 2);
  for (int i = 0; i < coordmat.nrow(); i++) {
    double x = coordmat(i,0);
    double y = coordmat(i,1);
    int col = round((x-xmin-(xres/2))/xres);
    int row = round((ymax-y-(yres/2))/yres);
    indexmat(i,0) = row;
    indexmat(i,1) = col;
  }
  return indexmat;
}

// [[Rcpp::export]]
NumericMatrix color_cpp(NumericMatrix rasterband, NumericMatrix coordvalmat, NumericVector extent, NumericVector res) {

  NumericMatrix newband = clone(rasterband);

  // Coordinates to index
  NumericMatrix coordmat = coordvalmat(_, Range(0,1));
  NumericMatrix indexmat = coord2index(coordmat, extent, res);

  // Color
  for (int i = 0; i < coordvalmat.nrow(); i++) {
    int row = indexmat(i, 0);
    int col = indexmat(i, 1);
    double value = coordvalmat(i,2);
    newband(row,col) = value;
  }

  return(newband);
}



