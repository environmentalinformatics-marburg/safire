### environmental stuff -----

## clear workspace 
rm(list = ls(all = TRUE))

## required packages
# devtools::install_github("MatMatt/MODIS", ref = "develop")
lib <- c("MODIS", "doParallel", "rworldmap")
jnk <- sapply(lib, function(x) library(x, character.only = TRUE))

## working directory
setwd("/media/fdetsch/modis_data")

## parallelization
cl <- makeCluster(detectCores() - 1)
registerDoParallel(cl)

## modis options
MODISoptions(localArcPath = "MODIS_ARC/", outDirPath = "MODIS_ARC/PROCESSED/", 
             outProj = "+init=epsg:4326")


### data download -----

## reference extent
ext <- subset(countriesCoarse, ADMIN.1 == "South Africa")

## download data in parallel
clc = getCollection("M*D14A1", forceCheck = TRUE)
tfs = runGdal("M*D14A1", extent = ext, job = "MCD14A1.006",
              collection = clc[[1]])


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
  
## deregister parallel backend
stopCluster(cl)
