# downloading and processing elevation(WC) and climate(CHELSA) data


# Load libraries ----------------------------------------------------------

rm(list=(ls())) # clear workspace
library(data.table)
library(terra)
library(sf)
library(geodata)





# Load data --------------------------------------------------------------------


atlas_years_cz = data.table(atlas = 5, 
                            years = c(2014:2017, 1985:1989),
                            period = c(2,2,2,2,1,1,1,1,1))
atlas_years_ny = data.table(atlas = 6, 
                            years = c(2000:2005,1980:1985),
                            period = c(2,2,2,2,2,2,1,1,1,1,1,1))
atlas_years_nz = data.table(atlas = 17, 
                            years = c(1999:2004,1969:1979),
                            period = c(2,2,2,2,2,2,1,1,1,1,1,1,1,1,1,1,1))
atlas_years_eu = data.table(atlas = 26, 
                            years = c(2013:2017, 1972:1995),
                            period = c(c(2,2,2,2,2), rep(1, length(1972:1995))))
atlas_years = rbind(atlas_years_cz, atlas_years_ny, atlas_years_nz, atlas_years_eu)





# Get geodata ---------------------------------------------------------------

## download once

geodata::elevation_global(path = "data/environment", res=0.5)



# we will get mean annula temperature, temperature seasonality, annual
# precipitation and precipitation seasonality. based on temp and pr.


years  <- unique(atlas_years[, years])
months <- sprintf("%02d", 1:12) 

base_url <- "https://os.unil.cloud.switch.ch/chelsa02/chelsa/global/monthly"

# check for downloaded already
# downloaded_clim = fread("downloaded_clim")
# years = years[!years %in% unique(unlist(regmatches(downloaded_clim$V1, gregexpr("[0-9]{4}", downloaded_clim$V1))))]

# chelsa monthly starts 
for(y in years[years>=1979]){
  for(m in months){
    
    ## Temperature
    tas_url <- sprintf(
      "%s/tas/%d/CHELSA_tas_%s_%d_V.2.1.tif",
      base_url, y, m, y)
    tas_dest <- sprintf(
      "CHELSA_tas_%s_%d_V.2.1.tif",
      m, y)
    
    ## Precipitation
    pr_url <- sprintf(
      "%s/pr/%d/CHELSA_pr_%s_%d_V.2.1.tif",
      base_url, y, m, y)
    pr_dest <- sprintf(
      "CHELSA_pr_%s_%d_V.2.1.tif",
      m, y)
    
    if(!file.exists(tas_dest)){
      try(
        download.file(tas_url, tas_dest, mode = "wb"),
        silent = TRUE
      )
    }
    if(!file.exists(pr_dest)){
      try(
        download.file(pr_url, pr_dest, mode = "wb"),
        silent = TRUE
      )
    }
  }
}



# missing years from daily ----
years <- 1969:1979
nz <- geodata::gadm("NZL", level = 0, path = tempdir())

base_url <- "https://os.unil.cloud.switch.ch/chelsa02/chelsa/global/daily/tas"

out_dir <- "CHELSA_annual"
tmp_dir <- "CHELSA_daily_tmp"

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tmp_dir, recursive = TRUE, showWarnings = FALSE)

terraOptions(
  tempdir = tmp_dir,
  threads = parallel::detectCores() - 1
)

for (y in years) {
  
  cat("\nProcessing", y, "\n")
  monthly_tas <- vector("list", 12)
  for (m in 1:12) {
    
    cat("  Month", m, "\n")
    
    ndays <- as.integer(format(
      seq(as.Date(sprintf("%d-%02d-01", y, m)),
          length = 2,
          by = "month")[2] - 1,
      "%d"))
    
    tas_list <- vector("list", ndays)
    
    #------------------------------------------------------
    # DOWNLOAD, CROP & MASK DAILY FILES
    #------------------------------------------------------
    
    for (d in 1:ndays) {
      
      date <- sprintf("%02d_%02d_%04d", d, m, y)
      
      file <- file.path(
        tmp_dir,
        sprintf("CHELSA_tas_%s_V.2.1.tif", date)
      )
      
      url <- sprintf(
        "%s/%d/CHELSA_tas_%s_V.2.1.tif",
        base_url, y, date
      )
      
      if (!file.exists(file)) {
        download.file(url, file, mode = "wb", quiet = TRUE)
      }
      
      r <- rast(file)
      r <- crop(r, nz)
      r <- mask(r, nz)
      
      tas_list[[d]] <- r
      
      file.remove(file)
    }
    
    #------------------------------------------------------
    # MONTHLY MEAN TEMPERATURE
    #------------------------------------------------------
    
    monthly_tas[[m]] <- app(rast(tas_list), mean)
    
    rm(tas_list)
    gc()
  }
  
  #--------------------------------------------------------
  # ANNUAL VARIABLES
  #--------------------------------------------------------
  
  tas_year <- rast(monthly_tas)
  
  # BIO1: Annual Mean Temperature (°C)
  BIO1 <- app(tas_year, mean) - 273.15
  names(BIO1) <- "BIO1"
  
  # BIO4: Temperature Seasonality
  BIO4 <- app(tas_year, sd) * 100
  names(BIO4) <- "BIO4"
  
  annual <- c(BIO1, BIO4)
  
  writeRaster(
    annual,
    filename = file.path(
      out_dir,
      sprintf("CHELSA_temperature_%d.tif", y)
    ),
    overwrite = TRUE
  )
  
  rm(monthly_tas, tas_year, BIO1, BIO4, annual)
  gc()
  
  cat("Finished", y, "\n")
}



# get annual variables
# set folders
terraOptions(tempdir=getwd())

# select atlas and years
#years  <- c(1985:1989, 2014:2017)
#cz <- geodata::gadm("CZ", path = "data/")
#dataset_id = 5

for(y in years){
  tas <- rast(sort(list.files(
    #"data/climate",
    pattern = paste0("_tas_", "[0-9]{2}", "_", y),
    full.names = TRUE)))
  #tas = crop(tas, cz)
  
  pr <- rast(sort(list.files(
    "data/climate",
    pattern = paste0("_pr_", "[0-9]{2}", "_", y),
    full.names = TRUE)))
  tap = crop(pr, cz)
  
  # MAT
  mat <- terra::app(tas, fun=mean)
  mat <- mat - 273.15 # to get from K to C
  MAT = round(mat,1)
  
  # TAP
  TAP <- sum(tap)
  
  # SEA_T
  SEA_T = terra::app(tas, fun=function(x){sd(x)*100})
  SEA_T = round(SEA_T, 1)
  
  # SEA_P
  SEA_P = terra::app(tap, fun=function(x){sd(x)/mean(x)})
  SEA_P = round(SEA_P, 2)
  
  writeRaster(MAT, filename=paste0("data/climate/annual/MAT_", y, "_", dataset_id, ".tif"),
              overwrite = TRUE,
              datatype = "FLT4S",
              gdal = c("COMPRESS=ZSTD","PREDICTOR=2"))
  writeRaster(TAP, filename=paste0("data/climate/annual/TAP_", y, "_", dataset_id,".tif"), 
              overwrite = TRUE,
              datatype = "FLT4S",
              gdal = c("COMPRESS=ZSTD","PREDICTOR=2"))
  writeRaster(SEA_T, filename=paste0("data/climate/annual/SEA_T_", y, "_", dataset_id,".tif"), 
              overwrite = TRUE,
              datatype = "FLT4S",
              gdal = c("COMPRESS=ZSTD","PREDICTOR=2"))
  writeRaster(SEA_P, filename=paste0("data/climate/annual/SEA_P_", y, "_", dataset_id,".tif"), 
              overwrite = TRUE,
              datatype = "FLT4S",
              gdal = c("COMPRESS=ZSTD","PREDICTOR=2"))
  
  terra::tmpFiles(remove=TRUE)
  
  cat(y, "\r")
}










# get annual variables
# set folders
terraOptions(tempdir=getwd())

# select atlas and years. avoid processing global maps, cropping is much faster

for(a in unique(atlas_years$atlas)){ # unique(atlas_years$atlas)
  years <- atlas_years[atlas == a, years]
  years <- years[years>=1979]
  
  avec <- terra::vect(paste0("all_scales_atlas_",a,".gpkg"))
  
  
  for(y in years){
    tas <- rast(sort(list.files(
      # "data/climate",
      pattern = paste0("_tas_", "[0-9]{2}", "_", y),
      full.names = TRUE)))
    tas = crop(tas, avec)
    
    tap <- rast(sort(list.files(
      # "data/climate",
      pattern = paste0("_pr_", "[0-9]{2}", "_", y),
      full.names = TRUE)))
    tap = crop(tap, avec)
    
    # MAT
    mat <- terra::app(tas, fun=mean)
    mat <- mat - 273.15 # to get from K to C
    MAT = round(mat,1)
    
    # TAP
    TAP <- sum(tap)
    
    # SEA_T
    SEA_T = terra::app(tas, fun=function(x){sd(x)*100})
    SEA_T = round(SEA_T, 1)
    
    # SEA_P
    SEA_P = terra::app(tap, fun=function(x){sd(x)/mean(x)})
    SEA_P = round(SEA_P, 2)
    
    writeRaster(MAT, filename=paste0(a,"_MAT_", y, ".tif"),
                overwrite = TRUE,
                datatype = "FLT4S",
                gdal = c("COMPRESS=ZSTD","PREDICTOR=2"))
    writeRaster(TAP, filename=paste0(a,"_TAP_", y, ".tif"), 
                overwrite = TRUE,
                datatype = "FLT4S",
                gdal = c("COMPRESS=ZSTD","PREDICTOR=2"))
    writeRaster(SEA_T, filename=paste0(a, "_SEA_T_", y, ".tif"), 
                overwrite = TRUE,
                datatype = "FLT4S",
                gdal = c("COMPRESS=ZSTD","PREDICTOR=2"))
    writeRaster(SEA_P, filename=paste0(a, "_SEA_P_", y, ".tif"), 
                overwrite = TRUE,
                datatype = "FLT4S",
                gdal = c("COMPRESS=ZSTD","PREDICTOR=2"))
    
    terra::tmpFiles(remove=TRUE)
    
    cat(y, "\r")
  }
  
}



# Process for difference between periods  --------------------------------------

atlas_years
vars = c("TAP", "MAT", "SEA_T", "SEA_P")

for(a in unique(atlas_years$atlas)){ # unique(atlas_years$atlas)
  years1 <- atlas_years[atlas == a & period==1, years]
  years1 <- years1[years1>=1979]
  years2 <- atlas_years[atlas == a & period==2, years]
  
  for(v in vars){
    t1 = dir(full.names = T, pattern = paste0("^", a, "_", v, "_", years1, ".tif", collapse = "|"))
    t2 = dir(full.names = T, pattern = paste0("^", a, "_", v, "_", years2, ".tif", collapse = "|"))
    
    t1 = rast(t1)
    t2 = rast(t2)
    
    # average per sampling period
    t1_mean = terra::app(t1, fun=mean)
    t2_mean = terra::app(t2, fun=mean)
    
    writeRaster(t1_mean, filename=paste0(a, "_",v,"_T1.tif"), 
                overwrite = TRUE,
                datatype = "FLT4S",
                gdal = c("COMPRESS=ZSTD","PREDICTOR=2")
    )
    writeRaster(t2_mean, filename=paste0(a, "_",v,"_T2.tif"), 
                overwrite = TRUE,
                datatype = "FLT4S",
                gdal = c("COMPRESS=ZSTD","PREDICTOR=2")
    )
    
    # difference
    dif = t2_mean - t1_mean
    names(dif)[which(names(dif)=="mean")] = v
    
    writeRaster(dif, filename=paste0(a, "_",v,"_diff.tif"),
                overwrite = TRUE,
                datatype = "FLT4S",
                gdal = c("COMPRESS=ZSTD","PREDICTOR=2")
    )
    terra::tmpFiles(remove=TRUE)
    cat(vars, "\r")
    
  }
}























