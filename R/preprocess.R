### environmental stuff -----

## clear workspace 
rm(list = ls(all = TRUE))

## required packages
# devtools::install_github("MatMatt/MODIS", ref = "develop")
lib <- c("MODIS", "doParallel", "rworldmap")
jnk <- sapply(lib, function(x) library(x, character.only = TRUE))

## working directory
setwd("/media/dogbert/modis_data")

## parallelization
cl <- makeCluster(detectCores() - 1)
registerDoParallel(cl)

## modis options
MODISoptions(localArcPath = "MODIS_ARC/", outDirPath = "MODIS_ARC/PROCESSED/")


### data download -----

# ## reference extent
# ext <- subset(countriesCoarse, ADMIN.1 == "South Africa")
# 
# ## download data in parallel
# foreach(product = c("MOD14A1", "MYD14A1"), .packages = "MODIS") %dopar%
#   runGdal(product, extent = ext, job = paste0(product, ".006"),
#           collection = getCollection(product, forceCheck = TRUE))


### data processing -----

## loop over products
lst_prd <- lapply(c("MOD14A1.006", "MYD14A1.006"), function(product) {
  
  dir_prd <- paste0("/media/dogbert/dev/data/", product)
  if (!dir.exists(dir_prd)) dir.create(dir_prd)
  
  ## crop images
  rst <- foreach(i = c("FireMask", "QA", "MaxFRP", "sample")) %do% {                                      
                       
    # list and import available files
    fls <- list.files(paste0(getOption("MODIS_outDirPath"), "/", product),
                      pattern = paste0(i, ".tif$"), full.names = TRUE)
    
    if (i %in% c("FireMask", "QA", "sample")) {
      rst <- raster::stack(fls)
    } else {
                       
      # apply scale factor
      dir_scl <- paste0(dir_prd, "/scl")
      if (!dir.exists(dir_scl)) dir.create(dir_scl)
      
      lst <- foreach(j = 1:length(fls), .packages = lib, 
                     .export = ls(envir = globalenv())) %dopar% {
        rst <- stack(fls[j])
        lst <- unstack(rst)
        rst <- stack(lapply(lst, function(k) {
          dataType(k) <- "INT4U"
          return(k)
        }))
        
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
                         format = "GTiff", overwrite = TRUE)
            }
          }
  
  rst_rcl <- stack(lst_rcl); rm(lst_rcl)
  
  
  ## create annual frequency images

  # indices
  dts <- extractDate(names(rst[[1]]))$inputLayerDates
  yrs <- substr(dts, 1, 4)
  ids <- as.numeric(as.factor(yrs))
  
  # target folder and files
  dir_yrs <- paste0(dir_prd, "/yrs")
  if (!dir.exists(dir_yrs)) dir.create(dir_yrs)
  
  fls_yrs <- sapply(unique(yrs), function(j) {
    gsub(dts[1], j, names(rst[[1]])[1])
  })
  fls_yrs <- paste0(dir_yrs, "/", fls_yrs, ".tif")
  
  for (j in seq(unique(ids))) {
    if (!file.exists(fls_yrs[j])) {
      rst_yr <- rst_rcl[[which(ids == unique(ids)[j])]]
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
