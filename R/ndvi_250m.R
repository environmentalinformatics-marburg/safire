### environment ----------------------------------------------------------------

## clear workspace
rm(list = ls(all = TRUE))

## load packages
lib <- c("parallel", "MODIS", "Rsenal")
Orcs::loadPkgs(lib)

## parallelization
cl <- makePSOCKcluster(detectCores() - 1)
jnk <- clusterEvalQ(cl, {library(MODIS); library(Rsenal)})

## 'MODIS' global settings
MODISoptions(localArcPath = "/media/fdetsch/data/MODIS_ARC/", 
             outDirPath = "/media/fdetsch/data/MODIS_ARC/PROCESSED/", 
             outProj = "+init=epsg:4326", 
             MODISserverOrder = c("LAADS", "LPDAAC"))


### preprocessing -----

## download .hdf files in parallel
library(rworldmap)
ref <- subset(countriesCoarse, ADMIN == "South Africa")
ref <- suppressWarnings(rgeos::gBuffer(ref, width = .1, quadsegs = 50L))

clusterExport(cl, "ref")
hdf <- parLapply(cl, c("MOD13A2", "MYD13A2"), function(product) {
  getHdf(product, extent = ref, collection = "006", 
         begin = "2008001", end = "2016366")
})

## extract required sds
tif <- vector("list", 2L); n <- 1
for (product in c("MOD13Q1", "MYD13Q1")) {
  tif[[n]] <- runGdal(product, extent = ref, collection = "006", 
                      job = paste0(product, ".006"), SDSstring = "101000000011", 
                      begin = "2008001", end = "2016366")
  n <- n + 1
}


### processing -----

## loop over products
setwd("/media/fdetsch/data/safire/")

lst_qc <- vector("list", 2L); n <- 1
for (product in c("MOD13A2", "MYD13A2")) {

  ## status message
  cat("Product '", product, "' is in, start processing.\n", sep = "")

  dir_prd <- paste0("data/", product, ".006")
  if (!dir.exists(dir_prd)) dir.create(dir_prd)

  ## import (and scale) images
  rst_sds <- lapply(c("NDVI", "pixel_reliability", "VI_Quality",
                      "composite_day_of_the_year"), function(i) {

    # list and import available files
    fls <- list.files(paste0(getOption("MODIS_outDirPath"), "/", product, ".006"),
                      pattern = paste0(i, ".tif$"), full.names = TRUE)
    rst <- raster::stack(fls)

    # scaling
    drs_scl <- paste0(dir_prd, "/scl")

    if (i %in% c("NDVI", "EVI")) {

      if (!dir.exists(drs_scl)) dir.create(drs_scl)
      fls_scl <- paste(getwd(), drs_scl, basename(fls), sep = "/")

      clusterExport(cl, c("fls_scl", "rst"), envir = environment())
      lst_scl <- parLapply(cl, 1:(raster::nlayers(rst)), function(j) {
        if (file.exists(fls_scl[j])) {
          raster::raster(fls_scl[j])
        } else {
          rst_scl <- rst[[j]] * 0.0001
          raster::writeRaster(rst_scl, filename = fls_scl[j],
                              format = "GTiff", overwrite = TRUE)
        }
      })

      return(raster::stack(lst_scl))

    } else {
      return(rst)
    }
  })


  ### quality control, step #1: -----
  ### discard clouds, snow/ice and filled pixels using 'pixel_reliability'

  ## status message
  cat("  Data import finished, starting with quality control #1.\n")

  dir_qc1 <- paste0(dir_prd, "/qc1")
  if (!dir.exists(dir_qc1)) dir.create(dir_qc1)

  ## perform quality check #1 for both NDVI and EVI
  fls_qc1 <- paste0(getwd(), "/", dir_qc1, "/", names(rst_sds[[1]]), ".tif")

  clusterExport(cl, c("fls_qc1", "rst_sds"), envir = environment())
  lst_qc1 <- parLapply(cl, 1:(raster::nlayers(rst_sds[[1]])), function(i) {
    if (file.exists(fls_qc1[i])) {
      raster(fls_qc1[i])
    } else {
      overlay(rst_sds[[1]][[i]], rst_sds[[2]][[i]], fun = function(x, y) {
        x[!y[] %in% c(0, 1)] <- NA
        return(x)
      }, filename = fls_qc1[i], overwrite = TRUE, format = "GTiff")
    }
  })

  rst_qc1 <- stack(lst_qc1); rm(lst_qc1)


  ### quality control, step #2: -----
  ### discard cloudy pixels based on 'state_250m' flags

  ## status message
  cat("  Quality control #1 finished, starting with quality control #2.\n")

  dir_qc2 <- paste0(dir_prd, "/qc2")
  if (!dir.exists(dir_qc2)) dir.create(dir_qc2)

  ## perform quality check #2 for both NDVI and EVI
  fls_qc2 <- paste0(getwd(), "/", dir_qc2, "/", names(rst_qc1), ".tif")

  clusterExport(cl, c("fls_qc2", "rst_qc1"), envir = environment())
  lst_qc2 <- parLapply(cl, 1:(raster::nlayers(rst_qc1)), function(i) {
    if (file.exists(fls_qc2[i])) {
      raster(fls_qc2[i])
    } else {
      overlay(rst_qc1[[i]], rst_sds[[3]][[i]], fun = function(x, y) {
        id <- sapply(y[], function(k) {
          bin <- number2binary(k, 16, TRUE)
          quality <- substr(bin, 15, 16)

          if (quality == "00") {
            return(TRUE)
          } else if (quality %in% c("10", "11")) {
            return(FALSE)
          } else {
            useful <- !substr(bin, 11, 14) %in% c("1101", "1110")
            aerosol <- substr(bin, 9, 10) != "11"
            adjacent <- substr(bin, 8, 8) == "0"
            mixed <- substr(bin, 6, 6) == "0"
            snow <- substr(bin, 2, 2) == "0"
            shadow <- substr(bin, 1, 1) == "0"

            all(useful, aerosol, adjacent, mixed, snow, shadow)
          }
        })

        x[!id] <- NA
        return(x)
      }, filename = fls_qc2[i], overwrite = TRUE, format = "GTiff")
    }
  })

  lst_qc[[n]] <- stack(lst_qc2); n <- n + 1
}


### whittaker smoother -----

## target folders and files
dir_prd <- "data/MCD13A2.006"
if (!dir.exists(dir_prd)) dir.create(dir_prd)

dir_wht <- paste0(dir_prd, "/whittaker")
if (!dir.exists(dir_wht)) dir.create(dir_wht)

## reorder layers
nms_qc2 <- do.call("c", lapply(lst_qc, names))
dts_qc2 <- extractDate(nms_qc2)$inputLayerDates
rst_qc2 <- stack(lst_qc)
rst_qc2 <- rst_qc2[[order(dts_qc2)]]
nms_qc2 <- nms_qc2[order(dts_qc2)]

## apply whittaker smoother
lst_wht <- whittaker.raster(rst_qc2, outDirPath = dir_wht,
                            overwrite = TRUE, format = "GTiff")

## write to disc
rst_wht <- stack(lst_wht)
names(rst_wht) <- gsub("MOD13A2", "MCD13A2", nms_qc2)
names(rst_wht) <- gsub("MYD13A2", "MCD13A2", names(rst_wht))
fls_wht <- paste0(getwd(), "/", dir_wht, "/", names(rst_wht), ".tif")

clusterExport(cl, c("rst_wht", "fls_wht"))
lst_wht <- parLapply(cl, 1:(raster::nlayers(rst_wht)), function(i) {
  rst <- rst_wht[[i]]
  rst[rst[] > 1 | rst[] < (-1)] <- NA
  
  writeRaster(rst, filename = fls_wht[i], format = "GTiff", overwrite = TRUE)
})

rst_wht <- stack(lst_wht)

## remove deprecated whittaker-related files
fls_old <- list.files(dir_wht, pattern = "yL5000", full.names = TRUE)
jnk <- file.remove(fls_old)

## deregister parallel backend
stopCluster(cl)
