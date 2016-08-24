#include <Rcpp.h>
#include <vector>

using namespace Rcpp;
using namespace std;

// [[Rcpp::export]]
NumericMatrix hittest_cpp(List rasterbands, NumericVector dim, NumericVector extent, NumericVector res, NumericVector polyX, NumericVector polyY, int polyCorners) {

  double xmin = extent[0];
  double ymax = extent[3];

  int nrow = dim[0];
  int ncol = dim[1];

  int nbands = rasterbands.size();

  double xres = res[0];
  double yres = res[1];

  vector<NumericMatrix> rbvec;
  for (int k = 0; k < nbands; k++) {
    NumericMatrix thisband = rasterbands[k];
    rbvec.push_back(thisband);
  }

  vector<double> constant;
  constant.reserve(polyCorners);
  vector<double> multiple;
  multiple.reserve(polyCorners);
  int i, j = polyCorners-1;
  bool oddNodes;
  double pxmin = polyX[0];
  double pxmax = polyX[0];
  double pymin = polyY[0];
  double pymax = polyY[0];
  double ty, tx;
  vector<double> hitx;
  vector<double> hity;
  vector<vector<double> > hitval;
  for (int k = 0; k < nbands; k++) {
    vector<double> val;
    hitval.push_back(val);
  }

  // Prepare polygon data
  for(i=0; i<polyCorners; i++) {

    // Polygon structure
    if(polyY[j]==polyY[i]) {
      constant[i]=polyX[i];
      multiple[i]=0;
    }
    else {
      constant[i]=polyX[i]-(polyY[i]*polyX[j])/(polyY[j]-polyY[i])+(polyY[i]*polyX[i])/(polyY[j]-polyY[i]);
      multiple[i]=(polyX[j]-polyX[i])/(polyY[j]-polyY[i]);
    }
    j=i;

    // Bounding box parameters
    if (polyX[i] < pxmin) {
      pxmin = polyX[i];
    }
    if (polyX[i] > pxmax) {
      pxmax = polyX[i];
    }
    if (polyY[i] < pymin) {
      pymin = polyY[i];
    }
    if (polyY[i] > pymax) {
      pymax = polyY[i];
    }
  }

  // Iterate through raster points and check whether in polygon
  int counter = 0;
  for (int trow=0; trow < nrow; trow++) {
    ty = ymax - trow*yres - (yres/2);
    for (int tcol=0; tcol < ncol; tcol++) {
      tx = xmin + tcol*xres + (xres/2);
      oddNodes=false;
      j=polyCorners-1;

      for (i=0; i<polyCorners; i++) {
        if (((polyY[i]< ty && polyY[j]>=ty)
               ||   (polyY[j]< ty && polyY[i]>=ty))) {
          oddNodes^=(ty*multiple[i]+constant[i]<tx);
        }
        j=i;
      }

      if (oddNodes) {
        hitx.push_back(tx);
        hity.push_back(ty);
        for (int k = 0; k < nbands; k++) {
          hitval[k].push_back(rbvec[k](trow, tcol));
        }
      }
      counter++;
    }
  }

  // Create output matrix
  NumericMatrix out (hitx.size(), 2+nbands);
  out(_,0) = (NumericVector)wrap(hitx);
  out(_,1) = (NumericVector)wrap(hity);
  for (int k = 0; k < nbands; k++) {
    out(_,2+k) = (NumericVector)wrap(hitval[k]);
  }
  return out;
}


// [[Rcpp::export]]
NumericMatrix unhit_cpp(NumericMatrix cmat, NumericVector polyX, NumericVector polyY, int polyCorners) {
  vector<double> constant;
  constant.reserve(polyCorners);
  vector<double> multiple;
  multiple.reserve(polyCorners);
  int i, j = polyCorners-1;
  bool oddNodes;
  double pxmin = polyX[0];
  double pxmax = polyX[0];
  double pymin = polyY[0];
  double pymax = polyY[0];
  double ty, tx;
  vector<double> hitx;
  vector<double> hity;
  vector<vector<double> > hitval;
  for (int k = 0; k < cmat.ncol()-2; k++) {
    vector<double> val;
    hitval.push_back(val);
  }

  // Prepare polygon data
  for(i=0; i<polyCorners; i++) {

    // Polygon structure
    if(polyY[j]==polyY[i]) {
      constant[i]=polyX[i];
      multiple[i]=0;
    }
    else {
      constant[i]=polyX[i]-(polyY[i]*polyX[j])/(polyY[j]-polyY[i])+(polyY[i]*polyX[i])/(polyY[j]-polyY[i]);
      multiple[i]=(polyX[j]-polyX[i])/(polyY[j]-polyY[i]);
    }
    j=i;

    // Bounding box parameters
    if (polyX[i] < pxmin) {
      pxmin = polyX[i];
    }
    if (polyX[i] > pxmax) {
      pxmax = polyX[i];
    }
    if (polyY[i] < pymin) {
      pymin = polyY[i];
    }
    if (polyY[i] > pymax) {
      pymax = polyY[i];
    }
  }

  // Iterate through raster points and check whether in polygon
  for (int p=0; p < cmat.nrow(); p++) {
    ty = cmat(p,1);
    tx = cmat(p,0);
    oddNodes=false;
    j=polyCorners-1;

    if (tx >= pxmin && tx <= pxmax && ty >= pymin && ty <= pymax) {
      for (i=0; i<polyCorners; i++) {
        if (((polyY[i]< ty && polyY[j]>=ty)
               ||   (polyY[j]< ty && polyY[i]>=ty))) {
          oddNodes^=(ty*multiple[i]+constant[i]<tx);
        }
        j=i;
      }
    }
    if (!oddNodes) {
      hitx.push_back(tx);
      hity.push_back(ty);
      for (int k = 0; k < cmat.ncol() - 2; k++) {
        hitval[k].push_back(cmat(p,2+k));
      }
    }
  }

  NumericMatrix out (hitx.size(), cmat.ncol());
  out(_,0) = (NumericVector)wrap(hitx);
  out(_,1) = (NumericVector)wrap(hity);
  for (int k = 0; k < cmat.ncol() - 2; k++) {
    out(_,2+k) = (NumericVector)wrap(hitval[k]);
  }

  return out;
}

