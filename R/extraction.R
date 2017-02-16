setwd("/media/fdetsch/data/safire/")

## ndvi data
fls_wht <- list.files("data/MCD13A2.006/whittaker", pattern = ".tif$", 
                      full.names = TRUE)

start <- grep("2010001", fls_wht); end <- grep("2014361", fls_wht)
fls_wht <- fls_wht[start:end]
rst_wht <- raster::stack(fls_wht)

## station data
spt <- readOGR("data/shp", "WeatherStations")
spt <- spTransform(spt, CRS = CRS(proj4string(ref)))
spt <- spt[ref, ]

## value extraction
val <- extract(rst_wht, spt, sp = TRUE)
val@data[, 3:ncol(val@data)] <- round(val@data[, 3:ncol(val@data)], 3)

saveRDS(val, "out/WeatherStationsNDVI.rds")
write.csv(val, "out/WeaterStationsNDVI.csv", row.names = FALSE, quote = FALSE)
