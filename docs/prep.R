
## Dependencies
library(velox)
library(sf)
library(raster)
library(rasterVis)

## Load data
nl.vx <- velox('/home/hunzikp/Data/viirs/SVDNB_npp_20170901-20170930_75N060W_vcmslcfg_v10_c201710041620/SVDNB_npp_20170901-20170930_75N060W_vcmslcfg_v10_c201710041620.avg_rade9.tif')
sui.sf <- st_read(dsn = '/home/hunzikp/Data/swisstopo/gemeinden/swissBOUNDARIES3D_LV03/SHAPEFILE_LV03_LN02', 
                  layer = 'swissBOUNDARIES3D_1_3_TLM_HOHEITSGEBIET')

## Get rid of lakes, etc
sui.sf <- sui.sf[!is.na(sui.sf$EINWOHNERZ),]

## Reproject
sui.sf <- st_transform(sui.sf, nl.vx$crs)

## Crop nighlights
nl.vx$crop(sui.sf)

## Save data
saveRDS(object = sui.sf, file = 'sui.Rda')
nl.vx$write(path = 'nl.tif', overwrite = TRUE)





## Extract using vx
vx.time <- system.time(ex.mat <- nl.vx$extract(sp = sui.sf, fun = function(x) sum(x, na.rm = TRUE)))

## Extract using raster
sui.sp <- as(st_zm(sui.sf), 'Spatial')
ras.time <- system.time(extract(x = nl.vx$as.RasterLayer(), sui.sp))

## Plot nightlights
nightPalette <- colorRampPalette(c("black", 
                                   rgb(4,6,84, maxColorValue=255),
                                   rgb(218,165,32, maxColorValue=255), 
                                   "white"), bias = 0.95)
nightTheme <- rasterTheme(region = nightPalette(20))
levelplot(nl.vx$as.RasterLayer(), 
          par.settings = nightTheme, margin=FALSE, xlab=NULL, ylab=NULL, 
          zscaleLog=TRUE, scales=list(draw=FALSE), colorkey = FALSE)

## Plot geometries
ar <- st_area(st_zm(sui.sf))
sui.sf$lnl <- log(as.vector(ex.mat)/as.vector(ar))
par(bg = grey(0.25), mar = rep(0, 4), xaxs = "i", yaxs = "i")
plot(sui.sf[,'lnl'], pal = nightPalette, main = NULL, key.pos = NULL)

