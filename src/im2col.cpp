#include <Rcpp.h>

using namespace Rcpp;


// [[Rcpp::export]]
NumericMatrix im2col_cpp(NumericMatrix rasterband, NumericVector dim,
                         int wrow, int wcol, int band, double padval, int rowframe, int colframe,
                         int rowstride, int colstride) {

  int nrow = dim[0];
  int ncol = dim[1];


  int nrowsteps = (int)ceil((double)nrow/(double)(wrow*rowstride));
  int ncolsteps = (int)ceil((double)ncol/(double)(wcol*colstride));

  NumericMatrix out(nrowsteps*ncolsteps, (wrow + 2*rowframe)*(wcol + 2*colframe));

  if (padval != 0) {
    int xsize = out.nrow() * out.ncol();
    for (int i = 0; i < xsize; i++) {
      out[i] = padval;
    }
  }

  int counter = 0;
  for (int i = 0; i < nrowsteps; i++) {
    for (int j = 0; j < ncolsteps; j++) {
      int iorigin = i*(wrow*rowstride) - rowframe;
      int jorigin = j*(wcol*colstride) - colframe;
      int wcounter = 0;
      for (int iw = iorigin; iw < iorigin+wrow+(2*rowframe); iw++) {
        for (int jw = jorigin; jw < jorigin+wcol+(2*colframe); jw++) {
          if (iw<nrow && iw >= 0 && jw<ncol && jw >= 0) {
            out(counter, wcounter) = rasterband(iw,jw);
          }
          wcounter++;
        }
      }
      counter++;
    }
  }

  return(out);
}


// [[Rcpp::export]]
NumericMatrix col2im_cpp(NumericMatrix rasterband, NumericVector dim,
                         NumericMatrix colmat, int wrow, int wcol, int band, int rowframe, int colframe,
                         int rowstride, int colstride) {

  NumericMatrix newband = clone(rasterband);

  int nrow = dim[0];
  int ncol = dim[1];

  int nrowsteps = (int)ceil((double)nrow/(double)(wrow*rowstride));
  int ncolsteps = (int)ceil((double)ncol/(double)(wcol*colstride));

  int counter = 0;
  for (int i = 0; i < nrowsteps; i++) {
    for (int j = 0; j < ncolsteps; j++) {
      int iorigin = i*(wrow*rowstride) - rowframe;
      int jorigin = j*(wcol*colstride) - colframe;
      int iorigin_noframe = i*(wrow*rowstride);
      int jorigin_noframe = j*(wcol*colstride);
      int imax_noframe = i*(wrow*rowstride) + wrow;
      int jmax_noframe = j*(wcol*colstride) + wcol;
      int wcounter = 0;
      for (int iw = iorigin; iw < iorigin+wrow+(2*rowframe); iw++) {
        for (int jw = jorigin; jw < jorigin+wcol+(2*colframe); jw++) {
          // Check whether iw-jw in image
          if (iw<nrow && iw >= 0 && jw<ncol && jw >= 0) {
            // Check whether iw-jw in patch
            if (iw >= iorigin_noframe && iw < imax_noframe && jw >= jorigin_noframe && jw < jmax_noframe){
              newband(iw,jw) = colmat(counter, wcounter);
            }
          }
          wcounter++;
        }
      }
      counter++;
    }
  }

  return(newband);
}
