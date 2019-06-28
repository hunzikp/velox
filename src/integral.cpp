#include <Rcpp.h>
#include <cmath>
using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix integralRaster(NumericMatrix rasterband) {

  int rows = rasterband.nrow();
  int cols = rasterband.ncol();



  for (int r = 0; r < rows; r++){
    for (int c = 0; c < cols; c++){
      if (c == 0){
        if (r != 0){
          rasterband(r, c) = rasterband(r, c) + rasterband(r - 1, c);
        }
      } else {
        if(r == 0){
          rasterband(r, c) = rasterband(r, c) + rasterband(r, c - 1);
        } else {
          rasterband(r, c) = rasterband(r, c) + rasterband(r - 1, c) +
            rasterband(r, c - 1) - rasterband(r - 1, c - 1);
        }
      }
    }
  }
  return  rasterband;
}

void fillFirstRows(int* firstRows, NumericMatrix kernel, int diameter){
  // group columns of kernel into boxes with same first and last row to make
  // even more use of integral image

  int prevFirstRow;
  int currValue;
  int prevValue;

  // iterate over cols to find first row of the kernel in each col;
  for (int c = 0; c < diameter; c++){

    // reset last value at beginning of each column
    prevValue = 0;

    // go through rows and check whether the value changed
    for (int r = 0; r < diameter; r++){

      // get the current value in the kernel matrix
      currValue = kernel(r,c);

      // if the value changed, then disc shape started/ended for this column
      if (prevValue != currValue){

        // if value changed from 0 to 1
        if (prevValue == 0){
          firstRows[c] = r;
          if (r != prevFirstRow) {
            prevFirstRow = r;
          }
          prevValue = currValue;
        }
      }
    }
    // return firstRows;
  }
}

void fillKernel(NumericMatrix kernel, int diameter){
  int radius = (diameter-1) / 2;
  int centX = radius;
  int centY = radius;
  int dX, dY;
  float dist;

  for (int y = 0; y < diameter; y++) {
    dY = abs(y - centY);
    for (int x = 0; x < diameter; x++){
      dX = abs(x - centX);
      dist = sqrt(pow(dX, 2) + pow(dY, 2));
      if(dist <= radius) {
        kernel(x,y) = 1;
      }
    }
  }
}


// [[Rcpp::export]]
NumericMatrix intFocalMean(NumericMatrix rasterband, int diameter) {

  // create output raster band
  int rows = rasterband.nrow();
  int cols = rasterband.ncol();
  NumericMatrix newband(rows, cols);

  // derive radius from diameter
  int radius;
  if (diameter % 2 == 1){
    radius = (diameter - 1) / 2;
  } else {
    radius = diameter / 2;
    diameter = diameter + 1;
  }

  // create disc shaped kernel
  NumericMatrix kernel(diameter, diameter);
  fillKernel(kernel, diameter);

  // find first rows in each column which are part of the disc kernel
  int firstRows[diameter];
  fillFirstRows(firstRows, kernel, diameter);

  for (auto x : firstRows) std::cout << x << " ";
  std::cout << std::endl;

  return kernel;

  // for (int r = 0; r < rows; r++){
  //   for (int c = 0; c < cols; c++){
  //
  //
  //
  //   }
  // }
  // return  rasterband;
}


std::list<int> getBoxCoordinates(int firstRows[]){
  std::list<int> ullr;

  return ullr;
}

/*** R
m <- matrix(c(0), nrow = 9, ncol = 9, byrow = TRUE)
m[5,5] <- 1

(i <- integralRaster(m))
(imn <- intFocalMean(i, 9))
raster::plot(raster::raster(imn))
# largeM <- matrix(runif(100000000), nrow = 10000)
# largeMint <- integralRaster(largeM)
*/
