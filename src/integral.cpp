#include <Rcpp.h>
#include <cmath>
#include <vector>
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


void fillKernel(NumericMatrix kernel, int diameter){
  // ancilliary function which uses the Pythagorean theorem to determine which
  // cell of the disc shaped focal kernel is within the disc shaped

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

void fillFirstRows(int* firstRows, NumericMatrix kernel, int diameter){
  // group columns of kernel into boxes with same first and last row to make
  // even more use of the integral image by creating rectangular boxes

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




std::vector<std::vector<int>> getBoxCoordinates(int * firstRows, int nRows){
  // anciliary function to get the upper left and lower right coordinates of
  // boxes. uses the first rows of the squared kernel which are occupied by the
  // circular shape if the filter to create groups of rectangular boxes which
  // start at the same line in the matrix. That is, the algorithm can make use
  // of the integral image and gain in efficiency.
  // This function returns a list of the upper left / lower right coordinates
  // of the individual boxes
  std::vector<std::vector<int>> ullr;

  int prevVal;          // previous value of firstRows, needed for detecting change
  int groupStartX = 0;  // remember the column in which the current group started
  std::vector<int> box;  // initialize vector for ullr coordinates of group


  for (int x = 0; x <= nRows; x++) {
    // at the first iteration set up group variable and remember the start of the
    // group in x direction
    if (x == 0){
      prevVal = firstRows[x];
      groupStartX = x;
    } else {
      // if the firstrow does not change (i.e. the kernel starts at the same
      // line as in the column before), just go on.

      if (firstRows[x] != prevVal){
      // if the kernel changes the row, create a new box for the group (i.e.
      // rectangle box) which has just finished.

        box.push_back(groupStartX); // left
        box.push_back(prevVal); // upper
        box.push_back(x - 1); // right
        box.push_back(nRows - 1 - prevVal); // lower

        // add the box to the the output
        ullr.push_back(box);
        box.clear();

        // reassign prevVal to currenct group's firstRow's
        prevVal = firstRows[x];
        // and save x position
        groupStartX = x;
      }
    }
  }

  return ullr;
}


// [[Rcpp::export]]
NumericMatrix intRadialFocalMean(NumericMatrix rasterband, int diameter) {
  // Implementation of a radial focal mean algorithm using itegral raster images.
  // This implementaiton is advantagous to classical focal operators, as it is
  // more efficient especially with large kernels.
  // With radial filters, one probably has to sacrifice some efficiency, though.

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

  // count number of pixels in kernel; needed for creating the mean
  int nPix = 0;
  for (int i = 0; i < diameter; i++)
    for (int j = 0; j < diameter; j++)
      if (kernel(i,j) == 1) nPix++;

  // find first rows in each column which are part of the disc kernel
  int firstRows[diameter];
  fillFirstRows(firstRows, kernel, diameter);

  // std::cout << "Firstrows: ";
  // for (auto x : firstRows) std::cout << x << " ";
  // std::cout << std::endl;

  // find rectangular boxes in the disc shaped kernel to speed up processing
  std::vector<std::vector<int>> ullr = getBoxCoordinates(firstRows, diameter);

  // convert ullr coordinates to values relative to the centering pixel of the kernel
  for (int x = 0; x < ullr.size(); x++) {
    for (int y = 0; y < ullr[x].size(); y++) {
      ullr[x][y] = ullr[x][y] - ((diameter - 1) / 2);
    }
  }


  // from here on copied and adapted from original velox/src/focal.cpp
  std::vector<double> v;
  int sum;

  int idim = (diameter-1)/2;
  int jdim = (diameter-1)/2;

  for (int i = 0; i < rows; i++) {
    int imin = i - idim;
    int imax = i + idim + 1;
    for (int j = 0; j < cols; j++) {
      int jmin = j - jdim;
      int jmax = j + jdim + 1;

      // for (int wi = imin; wi < imax; wi++) {
      //   for (int wj = jmin; wj < jmax; wj++) {
      //     if (wi >= 0 && wi < rows && wj >= 0 && wj < cols) {
      //       // v.push_back(rasterband(wi, wj));
      //     }
      //   }
      // }

      sum = accumulate(v.begin(),v.end(),0);


      newband(i,j) = sum / nPix;;
      v.clear();
    }
  }
  return kernel;
  // return newband;
}


/*** R
m <- matrix(c(0), nrow = 9, ncol = 9, byrow = TRUE)
m[5,5] <- 1

(i <- integralRaster(m))
(imn <- intRadialFocalMean(i, 9))
raster::plot(raster::raster(imn))
# largeM <- matrix(runif(100000000), nrow = 10000)
# largeMint <- integralRaster(largeM)
*/
