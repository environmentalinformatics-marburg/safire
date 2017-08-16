### environmental stuff -----

## clear workspace 
rm(list = ls(all = TRUE))

## required packages
# devtools::install_github("MatMatt/MODIS", ref = "develop")
lib <- c("MODIS", "doParallel", "rworldmap")
jnk <- sapply(lib, function(x) library(x, character.only = TRUE))

## working directory
repo = getwd()
setwd("/media/fdetsch/modis_data")

## parallelization
cl <- makeCluster(detectCores() - 1)
registerDoParallel(cl)

clusterExport(cl, "lib")
jnk = clusterEvalQ(cl, Orcs::loadPkgs(lib))

## options
options(stringsAsFactors = FALSE)

MODISoptions(localArcPath = "MODIS_ARC/", outDirPath = "MODIS_ARC/PROCESSED/", 
             outProj = "+init=epsg:4326")


### data download -----

## reference extent
ext <- subset(countriesCoarse, ADMIN.1 == "South Africa")

## download data
fun_tfs = function(file) {
  clc = MODIS::getCollection("M*D14A1", forceCheck = TRUE)
  tfs = MODIS::runGdal("M*D14A1", extent = ext, job = "MCD14A1.006",
                       collection = clc[[1]])
  saveRDS(tfs, file)
}

tfs = Orcs::ifMissing(file.path(repo, "inst/extdata/mcd14a1-tif.rds"), 
                      fun0 = readRDS, fun1 = fun_tfs, arg1 = "file")


### data processing -----

## loop over products
lst_prd <- lapply(names(tfs), function(product) {
  
  dir_prd <- paste0("/media/fdetsch/XChange/safire/", product)
  jnk = Orcs::ifMissing(dir_prd, file.path, dir.create, arg1 = "path")

  ## crop images
  rst <- foreach(i = 1:4) %do% {                                      
                       
    # list and import available files
    fls = sapply(tfs[[product]], "[[", i)

    if (length(grep("FireMask.tif$|QA.tif$|sample.tif$", fls[1])) == 1) {
      rst <- raster::stack(fls)
    } else {
                       
      # apply scale factor
      dir_scl <- paste0(dir_prd, "/scl")
      jnk = Orcs::ifMissing(dir_scl, file.path, dir.create, arg1 = "path")
      
      lst <- foreach(j = 1:length(fls), .packages = lib, 
                     .export = ls(envir = globalenv())) %dopar% {
        rst <- stack(fls[j])
        
        fls_scl <- paste0(dir_scl, "/", names(rst), ".tif")
        
        lst_scl <- lapply(1:nlayers(rst), function(k) {
          if (file.exists(fls_scl[k])) {
            raster::raster(fls_scl[[k]])
          } else {
            rst[[k]] <- rst[[k]] * 0.1
            raster::writeRaster(rst[[k]], filename = fls_scl[k], 
                                format = "GTiff", overwrite = TRUE)
          }
        })
        
        return(stack(lst_scl))
      }
      
      rst <- stack(lst)
    }
    
    return(rst)
  }

  ## reclassify images
  dir_rcl <- paste0(dir_prd, "/rcl")
  if (!dir.exists(dir_rcl)) dir.create(dir_rcl)
  
  fls_rcl <- paste0(dir_rcl, "/", names(rst[[1]]), ".tif")
  
  # reclassification matrix
  rcl <- matrix(c(0, 7, 0, 
                  7, 10, 1, 
                  10, 255, 0), ncol = 3, byrow = TRUE)
  
  lst_rcl <- foreach(j = 1:nlayers(rst[[1]]), .packages = lib, 
          .export = ls(envir = globalenv())) %dopar% {
            if (file.exists(fls_rcl[j])) {
              raster(fls_rcl[j])
            } else {
              reclassify(rst[[1]][[j]], rcl, include.lowest = TRUE, 
                         right = FALSE, filename = fls_rcl[j], 
                         overwrite = TRUE, datatype = "INT1U")
            }
          }
  
  rst_rcl <- stack(lst_rcl); rm(lst_rcl)
  
  
  ## create annual frequency images

  # indices
  nms = basename(sapply(tfs[[product]], "[[", 1))
  dts <- extractDate(nms)$inputLayerDates
  yrs <- substr(dts, 1, 4)

  # target folder and files
  dir_yrs <- paste0(dir_prd, "/yrs")
  if (!dir.exists(dir_yrs)) dir.create(dir_yrs)
  
  fls_yrs <- sapply(unique(yrs), function(j) {
    gsub(dts[1], j, nms[1])
  })
  fls_yrs <- paste0(dir_yrs, "/", fls_yrs, ".tif")
  
  for (j in seq(unique(yrs))) {
    if (!file.exists(fls_yrs[j])) {
      rst_yr <- rst_rcl[[which(yrs == unique(yrs)[j])]]
      calc(rst_yr, fun = function(x) {
        sum(x, na.rm = TRUE) / length(x)
      }, filename = fls_yrs[j], format = "GTiff", overwrite = TRUE)
      rm(rst_yr)
    }
  }
  
  rst_frq <- stack(fls_yrs)
  return(rst_frq)
})
  

### combined product -----

fun_hdf = function(file) {
  clc = MODIS::getCollection("M*D14A1", forceCheck = TRUE)
  fls = MODIS::getHdf("M*D14A1", extent = ext, collection = clc[[1]])
  saveRDS(fls, file)
  return(fls)
}

hdf = Orcs::ifMissing(file.path(repo, "inst/extdata/mcd14a1-hdf.rds"), 
                      fun0 = readRDS, fun1 = fun_hdf, arg1 = "file")

## loop over products
mcd_dts <- lapply(names(hdf), function(product) {
  
  ## dates from hdf files
  dpl = which(duplicated(dirname(hdf[[product]])))
  fls = hdf[[product]][-dpl]
  
  dates = do.call(c, parSapply(cl, fls, function(i) {
    nfo = suppressWarnings(rgdal::GDALinfo(i, returnScaleOffset = FALSE))
    mtd = attr(nfo, "mdata")
    dts = gsub("^Dates=", "", mtd[grep("^Dates=", mtd)])
    Orcs::unlistStrsplit(dts, " ")
  }))

  ## reclassified images
  rcl = list.files(paste0("/media/fdetsch/XChange/safire/", product, "/rcl"), 
                   pattern = ".tif$", full.names = TRUE)
  
  dat = data.frame(date = as.Date(dates), rcl, row.names = NULL)
  names(dat)[2] = product
  
  return(dat)
})

mcd_mrg = do.call(function(...) merge(..., by = "date", all = TRUE), mcd_dts)

## remove duplicated dates at turns of the year
dpl = which(duplicated(mcd_mrg$date))
mcd_mrg = mcd_mrg[-dpl, ]

mcd_fls = paste0("/media/fdetsch/XChange/safire/MCD14A1.006/MCD14A1.A", 
                 format(mcd_mrg$date, "%Y%j.FireMask.tif"))

clusterExport(cl, c("mcd_mrg", "mcd_fls"))
mcd_rst = parLapply(cl, 1:nrow(mcd_mrg), function(i) {
  if (is.na(mcd_mrg$MYD14A1.006[i])) {
    raster(mcd_mrg$MOD14A1.006[i])
  } else if (is.na(mcd_mrg$MOD14A1.006[i])) {
    raster(mcd_mrg$MYD14A1.006[i])
  } else {
    
    Orcs::ifMissing(mcd_fls[i], raster, function(filename) {
      mod = raster(mcd_mrg$MOD14A1.006[i])
      myd = raster(mcd_mrg$MYD14A1.006[i])
      overlay(mod, myd, fun = function(x, y) {
        x[y[] == 1] = 1L
        return(x)
      }, filename = mcd_fls[i], datatype = "INT1U")
    }, arg1 = "filename")
  }
})
mcd_rst = stack(mcd_rst)

## combined annual frequency images
dts <- extractDate(mcd_fls)$inputLayerDates
yrs <- substr(dts, 1, 4)

# target folder and files
dir_yrs <- paste0("/media/fdetsch/XChange/safire/MCD14A1.006/yrs")
if (!dir.exists(dir_yrs)) dir.create(dir_yrs)

fls_yrs <- sapply(unique(yrs), function(j) {
  gsub(dts[1], j, basename(mcd_fls[1]))
})
fls_yrs <- paste0(dir_yrs, "/", fls_yrs)

for (j in seq(unique(yrs))) {
  if (!file.exists(fls_yrs[j])) {
    rst_yr <- mcd_rst[[which(yrs == unique(yrs)[j])]]
    calc(rst_yr, fun = function(x) {
      sum(x, na.rm = TRUE) / length(x)
    }, filename = fls_yrs[j], format = "GTiff", overwrite = TRUE)
    rm(rst_yr)
  }
}

## deregister parallel backend
stopCluster(cl)
