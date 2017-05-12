####################################################
### NOTE
# For the PostGIS part of this script to run, you require
# an empty, PostGIS (>= 2.0) enabeld PostgreSQL DB (>=9.5) named
# 'gisdb' on your system, accessible via port 5432.
####################################################


####################################################
### INIT
####################################################
library(velox)
library(raster)
library(rgeos)
library(rbenchmark)
library(RPostgreSQL)
library(rgdal)

TEST_PGIS <- TRUE  # Switch off if no PostGIS DB available
PGIS_USER <- "hunzikp"  # Adjust accordingly
PGIS_PW <- "xxxxx"      # Adjust accordingly

####################################################
### Setup PGSQL connection
####################################################

if (TEST_PGIS) {
  drv <- dbDriver("PostgreSQL")
  con <- dbConnect(drv, dbname = "gisdb",
                   host = "localhost", port = 5432,
                   user = PGIS_USER, password = PGIS_PW)
}


####################################################
### Make testing data
####################################################

## Make large matrix
dim <- c(1000, 1000)
set.seed(0)
large.mat <- matrix(runif(prod(dim)), dim[1], dim[2])
extent <- c(0,1,0,1)
res <- 1/dim

## Make VeloxRaster and RasterLayer
large.vx <- velox(large.mat, extent, res, crs="")
large.rl <- raster(large.mat, xmn=extent[1], xmx=extent[2], ymn=extent[3], ymx=extent[4])

## Create SPDF
n <- 10
set.seed(0)
coords <- cbind(runif(n, extent[1], extent[2]), runif(n, extent[3], extent[4]))
sp <- SpatialPoints(coords)
spol <- gBuffer(sp, width=0.05, byid=TRUE)
spdf <- SpatialPolygonsDataFrame(spol, data.frame(id=1:length(spol)), FALSE)
writeOGR(spdf, "vignettes", "spdf", driver="ESRI Shapefile")

## Load data into PostGIS
if (TEST_PGIS) {
  ## Import raster into PostGIS and index
  dbSendQuery(con, "DROP TABLE IF EXISTS raster.large;")
  dbSendQuery(con, "DROP SCHEMA IF EXISTS raster;")
  dbSendQuery(con, "CREATE SCHEMA raster;")
  large.vx$write("vignettes/large.tif")
  system("raster2pgsql -s 4326 -I -C -M vignettes/large.tif -F -t 100x100 raster.large > importlarge.sql")
  system(paste0("psql -U ", PGIS_USER," -d gisdb -f importlarge.sql -h localhost -p 5432 -w"))
  file.remove("importlarge.sql")
  file.remove("vignettes/large.tif")

  ## Import SPDF into PostGIS
  dbSendQuery(con, "DROP TABLE IF EXISTS vector.spdf;")
  dbSendQuery(con, "DROP SCHEMA IF EXISTS vector;")
  dbSendQuery(con, "CREATE SCHEMA vector;")
  system("shp2pgsql -s 4326 -I  vignettes/spdf.shp vector.spdf > importspdf.sql")
  system(paste0("psql -U ", PGIS_USER, " -d gisdb -f importspdf.sql -h localhost -p 5432 -w"))
  file.remove("importspdf.sql")
  system("rm vignettes/spdf.*")
}


####################################################
### Prepare plot
####################################################

png("vignettes/benchmark.png", width=800, height=600)
layout(matrix(1:6, 2, 3, TRUE))


####################################################
### Benchmark extraction
####################################################

if (TEST_PGIS) {

  extract.postgis <- function() {
    query <- "SELECT
                id,
                SUM((ST_SUMMARYSTATS(ST_CLIP(rast, 1, geom, 0, FALSE))).sum) AS val
              FROM	raster.large,
              vector.spdf
              WHERE	ST_INTERSECTS(rast, geom)
              GROUP BY id
              ORDER BY id;"
    res <- dbSendQuery(con, query)
    ex.mat <- fetch(res, n = -1)
    dbClearResult(res)
  }

  ex.bm <- benchmark(velox=large.vx$extract(spdf, sum),
                     raster=extract(large.rl, spdf, fun=sum),
                     postgis=extract.postgis(),
                     replications=10,
                     columns=c("test", "replications", "elapsed", "relative"))

  ex.bm <- ex.bm[order(ex.bm$relative),]
  par(mar=c(5.1, 4.1, 7.1, 2.1))
  barplot(ex.bm$relative, names.arg=ex.bm$test,
          ylab="Relative Time Elapsed", cex.main=1.5)
  mtext("extract", line=5, cex=1.5, font=2)
  mtext(paste("velox ", round(ex.bm$relative[ex.bm$test=="raster"], 1), "x faster than raster", sep=""), line=3, cex=1.25)
  mtext(paste("velox ", round(ex.bm$relative[ex.bm$test=="postgis"], 1), "x faster than PostGIS", sep=""), line=1, cex=1.25)

} else {

  ex.bm <- benchmark(velox=large.vx$extract(spdf, sum),
                     raster=extract(large.rl, spdf, fun=sum),
                     replications=10,
                     columns=c("test", "replications", "elapsed", "relative"))

  ex.bm <- ex.bm[order(ex.bm$relative),]
  par(mar=c(5.1, 4.1, 7.1, 2.1))
  barplot(ex.bm$relative, names.arg=ex.bm$test,
          ylab="Relative Time Elapsed", cex.main=1.5)
  mtext("extract", line=5, cex=1.5, font=2)
  mtext(paste("velox ", round(ex.bm$relative[ex.bm$test=="raster"], 1), "x faster than raster", sep=""), line=3, cex=1.25)
}


####################################################
### Benchmark aggregation
####################################################


agg.velox <- function() {
  agg.vx <- large.vx$copy()
  agg.vx$aggregate(c(16,16), "mean")
}
ag.bm <- benchmark(velox=agg.velox(),
                   raster=aggregate(large.rl, fact=c(16,16), fun=mean, expand=FALSE),
                   replications=10,
                   columns=c("test", "replications", "elapsed", "relative"))

ag.bm <- ag.bm[order(ag.bm$relative),]
par(mar=c(5.1, 4.1, 6.1, 2.1))
barplot(ag.bm$relative, names.arg=ag.bm$test,
        ylab="Relative Time Elapsed", cex.main=1.5)
mtext("aggregate", line=4, cex=1.5, font=2)
mtext(paste("velox ", round(ag.bm$relative[ag.bm$test=="raster"], 1), "x faster than raster", sep=""), line=1.5, cex=1.25)


####################################################
### Benchmark cropping
####################################################


crop.sp <- spdf[1,]
crop.velox <- function() {
  crop.vx <- large.vx$copy()
  crop.vx$crop(crop.sp)
}
cr.bm <- benchmark(velox=crop.velox(),
                   raster=crop(large.rl, crop.sp),
                   replications=10,
                   columns=c("test", "replications", "elapsed", "relative"))

cr.bm <- cr.bm[order(cr.bm$relative),]
par(mar=c(5.1, 4.1, 6.1, 2.1))
barplot(cr.bm$relative, names.arg=cr.bm$test,
        ylab="Relative Time Elapsed", cex.main=1.5)
mtext("crop", line=4, cex=1.5, font=2)
mtext(paste("velox ", round(cr.bm$relative[cr.bm$test=="raster"], 1), "x faster than raster", sep=""), line=1.5, cex=1.25)


####################################################
### Benchmark median focal
####################################################

mfocal.velox <- function() {
  mfocal.vx <- large.vx$copy()
  mfocal.vx$medianFocal(5, 5)
}
weights <- matrix(1, 5, 5)

mf.bm <- benchmark(velox=mfocal.velox(),
                   raster=focal(large.rl, w=weights, fun=median, na.rm=TRUE, pad=TRUE, padValue=NA),
                   replications=10,
                   columns=c("test", "replications", "elapsed", "relative"))

mf.bm <- mf.bm[order(mf.bm$relative),]
par(mar=c(5.1, 4.1, 6.1, 2.1))
barplot(mf.bm$relative, names.arg=mf.bm$test,
        ylab="Relative Time Elapsed", cex.main=1.5)
mtext("median focal", line=4, cex=1.5, font=2)
mtext(paste("velox ", round(mf.bm$relative[mf.bm$test=="raster"], 1), "x faster than raster", sep=""), line=1.5, cex=1.25)




####################################################
### Benchmark sum focal
####################################################


weights <- matrix(1, 5, 5)
sfocal.velox <- function() {
  sfocal.vx <- large.vx$copy()
  sfocal.vx$sumFocal(weights, 1)
}

sf.bm <- benchmark(velox=sfocal.velox(),
                   raster=focal(large.rl, w=weights, fun=mean, na.rm=TRUE, pad=TRUE, padValue=NA),
                   replications=10,
                   columns=c("test", "replications", "elapsed", "relative"))

sf.bm <- sf.bm[order(sf.bm$relative),]
par(mar=c(5.1, 4.1, 6.1, 2.1))
barplot(sf.bm$relative, names.arg=sf.bm$test,
        ylab="Relative Time Elapsed", cex.main=1.5)
mtext("sum focal", line=4, cex=1.5, font=2)
mtext(paste("velox ", round(sf.bm$relative[sf.bm$test=="raster"], 1), "x faster than raster", sep=""), line=1.5, cex=1.25)



####################################################
### Benchmark rasterize
####################################################


rs.bm <- benchmark(velox=large.vx$rasterize(spdf, "id", 1, background=0),
                   raster=rasterize(spdf, large.rl, "id", background=0),
                   replications=10,
                   columns=c("test", "replications", "elapsed", "relative"))
rs.bm

rs.bm <- rs.bm[order(rs.bm$relative),]
par(mar=c(5.1, 4.1, 6.1, 2.1))
barplot(rs.bm$relative, names.arg=rs.bm$test,
        ylab="Relative Time Elapsed", cex.main=1.5)
mtext("rasterize", line=4, cex=1.5, font=2)
mtext(paste("velox ", round(rs.bm$relative[rs.bm$test=="raster"], 1), "x faster than raster", sep=""), line=1.5, cex=1.25)


dev.off()
