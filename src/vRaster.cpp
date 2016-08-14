#include <Rcpp.h>
#include <vector>
#include <numeric>
#include <math.h>
#include "gdal_priv.h"
#include "cpl_conv.h" 
#include "cpl_string.h"
#include "ogr_spatialref.h"

using namespace Rcpp;
using namespace std;

class vRaster {
  public:
    // Constructor
    vRaster();
    
    // Getter
    vector<NumericMatrix> getRasterBands();
    NumericMatrix getRasterBand(int band);
    NumericVector getResolution();
    NumericVector getExtent();
    NumericVector getDim();
    NumericMatrix getCoordinates();
    int getNBands();
    string getCRS();

    // Setter
    void setRasterBands(List x, NumericVector origin, NumericVector res, int nb, string incrs);
    
    // I/O
    void read(string path);
    void write(string path);
    
    // Crop Methods
    bool overlaps(NumericVector bbox);
    void crop(NumericVector bbox);
    vector<NumericMatrix> getCropBands(NumericVector cropbox);
    NumericVector getCropBox(NumericVector bbox);
    NumericVector getCropExtent(NumericVector cropbox);

    // Aggregate
    void aggregate(NumericVector factor, int aggtype);

    // Focal Functions
    void medianFocal(int wrow, int wcol, int band);
    void sumFocal(NumericMatrix weights, int wrow, int wcol, int band);
    void meanFocal(NumericMatrix weights, int wrow, int wcol, int band);
    
    // Vector Extraction
    NumericMatrix hittest(NumericVector polyX, NumericVector polyY, int polyCorners);
    NumericMatrix unhit(NumericMatrix cmat, NumericVector polyX, NumericVector polyY, int polyCorners);
    
  private:
    // Class Variables
    vector<NumericMatrix> rasterbands;
    int nbands;
    int nrow, ncol;
    double xmin, ymax;
    double xres, yres;
    string crs;
    
    // Math Helper Functions
    double median(vector<double> scores);
};


// Constructor
vRaster::vRaster() {
  nbands = 0;
  nrow = 0;
  ncol = 0;
  xmin = 0;
  ymax = 0;
  xres = 0;
  yres = 0;
  crs = "";
}


// Getters
vector<NumericMatrix> vRaster::getRasterBands() {
    return rasterbands;
}

NumericMatrix vRaster::getRasterBand(int band) {
    return rasterbands[band-1];
}

NumericVector vRaster::getResolution(){
  NumericVector res(2);
  res[0] = xres;
  res[1] = yres;
  return res;
}

NumericVector vRaster::getExtent(){
  NumericVector extent(4);
  double xmax = xmin + ncol*xres;
  double ymin = ymax - nrow*yres;
  extent[0] = xmin;
  extent[1] = xmax;
  extent[2] = ymin;
  extent[3] = ymax;
  return extent;
}

NumericVector vRaster::getDim(){
  NumericVector dim(2);
  dim[0] = nrow;
  dim[1] = ncol;
  return dim;
}

int vRaster::getNBands(){
  return nbands;
}

string vRaster::getCRS(){
  return crs;
}

NumericMatrix vRaster::getCoordinates(){
  NumericMatrix cmat = NumericMatrix(nrow*ncol, 2);
  int counter = 0;
  for (int i=0; i<nrow; i++) {
    for (int j=0; j<ncol; j++) {
      cmat(counter, 0) =  xmin + j*xres + (xres/2);
      cmat(counter, 1) =  ymax - i*yres - (yres/2);
      counter++;
    }
  }
  return cmat;
}


// Setters
void vRaster::setRasterBands(List x, NumericVector origin, NumericVector res, int nb, string incrs) {
    rasterbands.clear();
    for (int k=0; k<nb; k++) {
      NumericMatrix thisband = x[k];
      rasterbands.push_back(thisband);
    }
    nrow = rasterbands[0].nrow();
    ncol = rasterbands[0].ncol();
    xmin = origin[0];
    ymax = origin[1];
    xres = res[0];
    yres = res[1];
    nbands = nb;
    crs = incrs;
}


// Input/Output
void vRaster::read(string path) {
    GDALDataset  *poDataset;
    GDALAllRegister();
    poDataset = (GDALDataset *) GDALOpen( path.c_str(), GA_ReadOnly );
    
    // Get raster specs
    double adfGeoTransform[6];
    poDataset->GetGeoTransform( adfGeoTransform );
    xmin = adfGeoTransform[0];
    xres = adfGeoTransform[1];
    ymax = adfGeoTransform[3];
    yres = -1*adfGeoTransform[5];
    
    // Get Proj4
    OGRSpatialReference oSRS;
    const char *pszSRS_WKT = NULL;
    pszSRS_WKT = poDataset->GetProjectionRef();
    char * pszSRS_WKTchar = strdup(pszSRS_WKT);
    oSRS.importFromWkt(&pszSRS_WKTchar);
    char *pszSRS_Proj4 = NULL;
    oSRS.exportToProj4(&pszSRS_Proj4);
    crs = string(pszSRS_Proj4);

    // Read bands
    nbands = poDataset->GetRasterCount();
    
    GDALRasterBand  *poBand;
    poBand = poDataset->GetRasterBand( 1 );
    int   nXSize = poBand->GetXSize();
    int   nYSize = poBand->GetYSize();
    ncol = nXSize;
    nrow = nYSize;
    
    for (int k = 1; k<=nbands; k++) {
      poBand = poDataset->GetRasterBand( k );
      float *pafScanline;
      NumericMatrix newMatrix(Dimension(nYSize, nXSize));
      pafScanline = (float *) CPLMalloc(sizeof(float)*nXSize*nYSize);
      poBand->RasterIO( GF_Read, 0, 0, nXSize, nYSize, 
                    pafScanline, nXSize, nYSize, GDT_Float32 , 
                    0, 0 );
  
      int counter = 0;
      for (int i = 0; i < nYSize; i++) {
        for (int j = 0; j < nXSize; j++) {
          newMatrix(i, j) = static_cast<double>(pafScanline[counter]);
          counter++;
        }
      }
      rasterbands.push_back(newMatrix);
      CPLFree(pafScanline);
    }
    
    GDALClose((GDALDatasetH)poDataset);
}

void vRaster::write(string path) {
  // Create file
  const char *pszFormat = "GTiff";
  GDALDriver *poDriver;
  char **papszMetadata;
  poDriver = GetGDALDriverManager()->GetDriverByName(pszFormat);
  GDALDataset *poDstDS;       
  char **papszOptions = NULL;
  poDstDS = poDriver->Create(path.c_str(), ncol, nrow, nbands, GDT_Float32, 
                              papszOptions);
  double adfGeoTransform[6] = { xmin, xres, 0, ymax, 0, -yres};
  OGRSpatialReference oSRS;
  char *pszSRS_WKT = NULL;
  GDALRasterBand *poBand;
  poDstDS->SetGeoTransform( adfGeoTransform );

  // Set CRS
  oSRS.importFromProj4(crs.c_str());
  oSRS.exportToWkt( &pszSRS_WKT );
  poDstDS->SetProjection(pszSRS_WKT);
  CPLFree(pszSRS_WKT);
  
  // Write each band to file
  for (int k = 1; k <= nbands; k++) {
    float* rasterData = new float[ncol*nrow];
    NumericMatrix rasterband = rasterbands[k-1];
    int counter = 0;
    for (int i = 0; i < nrow; i++) {
      for (int j = 0; j < ncol; j++) {
        rasterData[counter] = rasterband(i, j);
        counter++;
      }
    }
  
    poBand = poDstDS->GetRasterBand(k);
    poBand->RasterIO( GF_Write, 0, 0, ncol, nrow, 
                    rasterData, ncol, nrow, GDT_Float32, 0, 0 );    
  }
  
  // Clean up
  GDALClose( (GDALDatasetH) poDstDS );
}


// Crop
bool vRaster::overlaps(NumericVector bbox) {
  double xmax = xmin + ncol*xres;
  double ymin = ymax - nrow*yres;
  double xminbb = bbox[0];
  double xmaxbb = bbox[1];
  double yminbb = bbox[2];
  double ymaxbb = bbox[3];
  
  if (xmaxbb > xmin && xminbb < xmax && ymaxbb > ymin && yminbb < ymax) {
    return true;
  } else {
    return false;
  }
}

void vRaster::crop(NumericVector bbox){
  NumericVector cropbox = getCropBox(bbox);
  NumericVector cropextent = getCropExtent(cropbox);
  rasterbands = getCropBands(cropbox);
  
  nrow = rasterbands[0].nrow();
  ncol = rasterbands[0].ncol();
  xmin = cropextent[0];
  ymax = cropextent[3];
}

NumericVector vRaster::getCropBox(NumericVector bbox) {
  double minrow, maxrow, mincol, maxcol;
  
  double xmax = xmin + ncol*xres;
  double ymin = ymax - nrow*yres;
  
  double xminbb = bbox[0];
  double xmaxbb = bbox[1];
  double yminbb = bbox[2];
  double ymaxbb = bbox[3];
  
  // Get crop dimensions
  if (xminbb > xmin) {
    mincol = floor((xminbb - xmin)/xres);
  } else {
    mincol = 0;
  }
  if (xmaxbb < xmax) {
    maxcol = ceil((xmaxbb - xmin)/xres);
  } else {
    maxcol = ncol;
  }
  if (yminbb > ymin) {
    maxrow = ceil((ymax - yminbb)/yres);
  } else {
    maxrow = nrow;
  }
  if (ymaxbb < ymax) {
    minrow = floor((ymax - ymaxbb)/yres);
  } else {
    minrow = 0;
  }
  // Adjust maxrow/maxcol for counting from zero
  maxrow = maxrow - 1;
  maxcol = maxcol - 1;
  // Output
  NumericVector cropbox = NumericVector::create(minrow, maxrow, mincol, maxcol);
  return cropbox;
}

vector<NumericMatrix> vRaster::getCropBands(NumericVector cropbox) {
  vector<NumericMatrix> cropbands;
  for (int k=0; k < nbands; k++) {
    NumericMatrix thisband = rasterbands[k];
    cropbands.push_back(thisband( Range(cropbox[0], cropbox[1]), Range(cropbox[2], cropbox[3])));
  } 
  return cropbands;
}


NumericVector vRaster::getCropExtent(NumericVector cropbox) {
  double xmax = xmin + ncol*xres;
  double ymin = ymax - nrow*yres;
  
  double xminn = xmin + cropbox[2]*xres;
  double xmaxn = xmin + cropbox[3]*xres;
  double yminn = ymax - cropbox[1]*yres;
  double ymaxn = ymax - cropbox[0]*yres;
  
  NumericVector cropextent =NumericVector::create(xminn, xmaxn, yminn, ymaxn);
  return cropextent;
}


// Aggregate
void vRaster::aggregate(NumericVector factor, int aggtype) {
  // aggtype 0: sum <<- default if aggtype outside allowed range
  // aggtype 1: mean
  // aggtype 2: min
  // aggtype 3: max
  // aggtype 4: median
  int nrown = (int)floor(nrow/factor[1]);
  int ncoln = (int)floor(ncol/factor[0]);
  int startrow, endrow, startcol, endcol;
  vector<NumericMatrix> newrasterbands;
  
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
    newrasterbands.push_back(newband);
  }
  
  xres = xres*factor[0];
  yres = yres*factor[1];
  nrow = nrown;
  ncol = ncoln;
  rasterbands = newrasterbands;
}


// Focal Methods
void vRaster::medianFocal(int wrow, int wcol, int band) {
  int i, j, mini, maxi, minj, maxj;
  int nbsize;
  int k;
  vector<double> v;
  NumericMatrix thisband = rasterbands[band-1];
  NumericMatrix newband(nrow, ncol);
  
  int idim = (wrow-1)/2;
  int jdim = (wcol-1)/2;
  
  for (i = 0; i < nrow; i++) {
    if (i - idim < 0) {
      mini = 0;
    } else {
      mini = i-idim;
    }
    if (i + idim >= nrow) {
      maxi = nrow-1;
    } else {
      maxi = i+idim;
    }
    for (j=0; j < ncol; j++) {
      if (j - jdim < 0) {
        minj = 0;
      } else {
        minj = j-jdim;
      }
      if (j + jdim >= ncol) {
        maxj = ncol-1;
      } else {
        maxj = j+jdim;
      }
      
      NumericMatrix nb = thisband(Range(mini,maxi), Range(minj,maxj));
      nbsize = (maxi-mini+1)*(maxj-minj+1);
      for (k = 0; k < nbsize; k++) {
        v.push_back(nb[k]);
      }
      newband(i,j) = median(v);
      v.clear();
    }
  }
  rasterbands[band-1] = newband;
}

void vRaster::sumFocal(NumericMatrix weights, int wrow, int wcol, int band) {
  int i, j, mini, maxi, minj, maxj, ki, kj, wi, wj;
  int nbsize;
  int k;
  double cusum;
  int counter;
  double w;
  NumericMatrix thisband = rasterbands[band-1];
  NumericMatrix newband(nrow, ncol);
  
  int idim = (wrow-1)/2;
  int jdim = (wcol-1)/2;
  
  for (i = 0; i < nrow; i++) {
    if (i - idim < 0) {
      mini = 0;
    } else {
      mini = i-idim;
    }
    if (i + idim >= nrow) {
      maxi = nrow-1;
    } else {
      maxi = i+idim;
    }
    for (j=0; j < ncol; j++) {
      if (j - jdim < 0) {
        minj = 0;
      } else {
        minj = j-jdim;
      }
      if (j + jdim >= ncol) {
        maxj = ncol-1;
      } else {
        maxj = j+jdim;
      }
      
      cusum = 0;
      for (ki = mini; ki < maxi; ki++) {
        for (kj = minj; kj < maxj; kj++) {
          wi = ((wrow-1)/2) + (ki - i);
          wj = ((wcol-1)/2) + (kj - j);
          w = weights(wi, wj);
          cusum = cusum + thisband(ki, kj)*w;
        }
      }
      newband(i,j) = cusum;
    }
  }
  rasterbands[band-1] = newband;
}

void vRaster::meanFocal(NumericMatrix weights, int wrow, int wcol, int band) {
  int i, j, mini, maxi, minj, maxj, ki, kj, wi, wj;
  int nbsize;
  int k;
  double cusum;
  double wsum;
  int counter;
  double w;
  NumericMatrix thisband = rasterbands[band-1];
  NumericMatrix newband(nrow, ncol);

  int idim = (wrow-1)/2;
  int jdim = (wcol-1)/2;
  
  for (i = 0; i < nrow; i++) {
    if (i - idim < 0) {
      mini = 0;
    } else {
      mini = i-idim;
    }
    if (i + idim >= nrow) {
      maxi = nrow-1;
    } else {
      maxi = i+idim;
    }
    for (j=0; j < ncol; j++) {
      if (j - jdim < 0) {
        minj = 0;
      } else {
        minj = j-jdim;
      }
      if (j + jdim >= ncol) {
        maxj = ncol-1;
      } else {
        maxj = j+jdim;
      }
      
      cusum = 0;
      wsum = 0;
      for (ki = mini; ki < maxi; ki++) {
        for (kj = minj; kj < maxj; kj++) {
          wi = ((wrow-1)/2) + (ki - i);
          wj = ((wcol-1)/2) + (kj - j);
          w = weights(wi, wj);
          cusum = cusum + thisband(ki, kj)*w;
          wsum = wsum + w;
        }
      }
      newband(i,j) = cusum/wsum;
    }
  }
  rasterbands[band-1] = newband;
}


// Vector Extraction
NumericMatrix vRaster::hittest(NumericVector polyX, NumericVector polyY, int polyCorners) {
  double constant[polyCorners];
  double multiple[polyCorners];
  int   i, p, j = polyCorners-1;
  bool  oddNodes;
  double pxmin = polyX[0];
  double pxmax = polyX[0];
  double pymin = polyY[0];
  double pymax = polyY[0];
  double ty, tx;
  double xminn, xmaxn, yminn, ymaxn;
  vector<NumericMatrix> cropbands;
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
  
  // Crop raster (if there is an overlap)
  NumericVector bbox = NumericVector::create(pxmin, pxmax, pymin, pymax);
  bool op = overlaps(bbox);
  if (op == true) {
    NumericVector cropbox = getCropBox(bbox);
    cropbands = getCropBands(cropbox);
    NumericVector cropextent = getCropExtent(cropbox);
    xminn = cropextent[0];
    xmaxn = cropextent[1];
    yminn = cropextent[2];
    ymaxn = cropextent[3];
  } else {
    NumericMatrix out(0,3);
    return out;
  }
  
  // Iterate through cropped raster points and check whether in polygon
  int counter = 0;
  for (int trow=0; trow < cropbands[0].nrow(); trow++) {
    ty = ymaxn - trow*yres - (yres/2);
    for (int tcol=0; tcol < cropbands[0].ncol(); tcol++) {
      tx = xminn + tcol*xres + (xres/2);
      oddNodes=false;
    	j=polyCorners-1;
      
      for (i=0; i<polyCorners; i++) {
    		if ((polyY[i]< ty && polyY[j]>=ty
    		||   polyY[j]< ty && polyY[i]>=ty)) {
      		oddNodes^=(ty*multiple[i]+constant[i]<tx); 
        }
        j=i; 
      }
      
      if (oddNodes) {
        hitx.push_back(tx);
        hity.push_back(ty);
        for (int k = 0; k < nbands; k++) {
          hitval[k].push_back(cropbands[k](trow, tcol));
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

NumericMatrix vRaster::unhit(NumericMatrix cmat, NumericVector polyX, NumericVector polyY, int polyCorners) {
  double constant[polyCorners];
  double multiple[polyCorners];
  int   i, p, j= polyCorners-1;
  bool  oddNodes;
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
  
  // Iterate through cropped raster points and check whether in polygon
  for (int p=0; p < cmat.nrow(); p++) {
    ty = cmat(p,1);
    tx = cmat(p,0);
    oddNodes=false;
    j=polyCorners-1;
    
    if (tx >= pxmin && tx <= pxmax && ty >= pymin && ty <= pymax) {
      for (i=0; i<polyCorners; i++) {
    		if ((polyY[i]< ty && polyY[j]>=ty
    		||   polyY[j]< ty && polyY[i]>=ty)) {
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


// Math Helper Functions
double vRaster::median(vector<double> scores) {
  double median;
  size_t size = scores.size();

  sort(scores.begin(), scores.end());

  if (size  % 2 == 0) {
      median = (scores[size / 2 - 1] + scores[size / 2]) / 2;
  }
  else {
      median = scores[size / 2];
  }

  return median;
}


// Rcpp Module Exposure
RCPP_MODULE(vRaster_module) {
  class_<vRaster>( "vRaster" )
  .constructor()
  .method( "getRasterBands", &vRaster::getRasterBands)
  .method( "getRasterBand", &vRaster::getRasterBand)
  .method( "getResolution", &vRaster::getResolution)
  .method( "getExtent", &vRaster::getExtent)
  .method( "getDim", &vRaster::getDim)
  .method( "getNBands", &vRaster::getNBands)
  .method( "getCRS", &vRaster::getCRS)
  .method( "getCropBands", &vRaster::getCropBands)
  .method( "getCropBox", &vRaster::getCropBox)
  .method( "getCropExtent", &vRaster::getCropExtent)
  .method( "getCoordinates", &vRaster::getCoordinates)
  .method( "overlaps", &vRaster::overlaps)
  .method( "setRasterBands", &vRaster::setRasterBands)
  .method( "read", &vRaster::read)
  .method( "write", &vRaster::write)
  .method( "crop", &vRaster::crop)
  .method( "aggregate", &vRaster::aggregate)
  .method( "medianFocal", &vRaster::medianFocal)
  .method( "sumFocal", &vRaster::sumFocal)
  .method( "meanFocal", &vRaster::meanFocal)
  .method( "hittest", &vRaster::hittest)
  .method( "unhit", &vRaster::unhit)
  ;
}
