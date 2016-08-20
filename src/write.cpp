#include <Rcpp.h>
#include <math.h>
#include "gdal_priv.h"
#include "cpl_conv.h"
#include "cpl_string.h"
#include "ogr_spatialref.h"

using namespace Rcpp;
using namespace std;

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
