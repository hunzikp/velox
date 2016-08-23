#include <Rcpp.h>
#include <math.h>
#include "gdal_priv.h"
#include "cpl_conv.h"
#include "cpl_string.h"
#include "ogr_spatialref.h"

using namespace Rcpp;
using namespace std;

// [[Rcpp::export]]
List readVelox(std::string path) {

  Rcpp::List out(6);

  GDALDataset  *poDataset;
  GDALAllRegister();
  poDataset = (GDALDataset *) GDALOpen( path.c_str(), GA_ReadOnly );

  // Get raster specs
  double adfGeoTransform[6];
  poDataset->GetGeoTransform( adfGeoTransform );
  double xmin = adfGeoTransform[0];
  double xres = adfGeoTransform[1];
  double ymax = adfGeoTransform[3];
  double yres = -1*adfGeoTransform[5];

  NumericVector res = NumericVector::create(xres, yres);
  out[0] = res;
  NumericVector origin = NumericVector::create(xmin, ymax);
  out[1] = origin;

  // Get Proj4
  OGRSpatialReference oSRS;
  const char *pszSRS_WKT = NULL;
  pszSRS_WKT = poDataset->GetProjectionRef();
  char * pszSRS_WKTchar = strdup(pszSRS_WKT);
  oSRS.importFromWkt(&pszSRS_WKTchar);
  char *pszSRS_Proj4 = NULL;
  oSRS.exportToProj4(&pszSRS_Proj4);
  string crs = string(pszSRS_Proj4);

  out[2] = crs;

  // Read bands
  int nbands = poDataset->GetRasterCount();
  out[3] = nbands;

  GDALRasterBand  *poBand;
  poBand = poDataset->GetRasterBand( 1 );
  int   nXSize = poBand->GetXSize();
  int   nYSize = poBand->GetYSize();
  int ncol = nXSize;
  int nrow = nYSize;

  NumericVector dim = NumericVector::create(nrow, ncol);
  out[4] = dim;

  Rcpp::List rasterbands(nbands);
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
    rasterbands[k-1] = newMatrix;
    CPLFree(pafScanline);
  }
  out[5] = rasterbands;
  GDALClose((GDALDatasetH)poDataset);

  return(out);
}


// [[Rcpp::export]]
void writeVelox(std::string path, List rasterbands, NumericVector dim, NumericVector extent, NumericVector res, std::string crs) {

  int nrow = dim[0];
  int ncol = dim[1];
  double xmin = extent[0];
  double ymax = extent[3];
  double xres = res[0];
  double yres = res[1];
  int nbands = rasterbands.size();

  // Create file
  GDALAllRegister();
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
