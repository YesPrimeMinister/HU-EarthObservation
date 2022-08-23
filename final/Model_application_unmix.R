# =============================================================================
# GCG - EARTH OBSERVATION SUMMER 2022
# Grassland droughts
# =============================================================================
# load required packages
l_packages <- c("assertthat","GGally", "ggplot2", "tidyr", "dplyr","raster", "e1071", 
                "rgdal", "viridisLite", "randomForest", "ggpubr", "Metrics", 'terra')
for (p in l_packages){
  if(! p %in% installed.packages()){
    install.packages(p, dependencies = TRUE)
  }
  library(p, character.only = T)
}

# =============================================================================
# Set working directory
dir_main <- 'data'
setwd(dir_main)

# =============================================================================
# Load SVR models
load(file='models/svr_pv_500_samples.txt')
load(file='models/svr_npv_500_samples.txt')
load(file='models/svr_soil_500_samples.txt')

# =============================================================================
# Define function to apply models and export the result to out_dir
unmix <- function(raster_list, out_dir, model_pv, model_npv, model_soil){
  
  "Function to unmix all sentinel-2 observations from FORCE and export rasters of individual observation dates.
  
  raster_list:        (list) Input filenames of Sentinel-2 observations exported from FORCE (e.g.: individual rasters represent spectral bands) 
  out_dir:            (char) Directory path where to export the results.
  model_pv:           (list) Model for predicting the fraction of photosynthetic vegetation.
  model_npv:          (list) Model for predicting the fraction of non-photosynthetic vegetation.
  model_soil:         (list) Model for predicting the fraction of soils.
  "
  
  # number of observations = number of bands in the input rasters
  n_bands <- nlayers(stack(raster_list[[1]]))
  # dates of sentinel-2 observations
  obs_dates <- names(stack(raster_list[[1]]))
  # spectral band names derived from the trained model
  band_names <- dimnames(model_pv$SV)[[2]]
  
  print(paste('Processing observation # out of ', n_bands, ':', sep=''))
  
  # unmix at every observation in time
  for (obs in 1:n_bands) {
    print(paste(obs, '/', n_bands, sep=''))
    
    # read inidvidual spectral bands from files
    blue  <- raster(raster_list[1], band = obs)
    green <- raster(raster_list[3], band = obs)
    red   <- raster(raster_list[8], band = obs)
    re1   <- raster(raster_list[5], band = obs)
    re2   <- raster(raster_list[6], band = obs)
    re3   <- raster(raster_list[7], band = obs)
    bnir  <- raster(raster_list[2], band = obs)
    nir   <- raster(raster_list[4], band = obs)
    swir1 <- raster(raster_list[9], band = obs)
    swir2 <- raster(raster_list[10], band = obs)
    
    # stack spectral bands into a monotemporal, multispectral raster stack
    stack_monotemp <- stack(blue, green, red, re1, re2, re3, bnir, nir, swir1, swir2)
    stack_monotemp <- stack_monotemp / 10000
    names(stack_monotemp) <- band_names
    
    # apply unmixing models to the monotemporal raster
    pred_pv   <- predict(stack_monotemp, model_pv,   type='response', na.omit=T)
    pred_npv  <- predict(stack_monotemp, model_npv,  type='response', na.omit=T)
    pred_soil <- predict(stack_monotemp, model_soil, type='response', na.omit=T)
    
    # add raster band name
    names(pred_pv)   <- 'pv'
    names(pred_npv)  <- 'npv'
    names(pred_soil) <- 'soil'
    
    # remove values above 1 / below 0
    pred_pv   <- clamp(pred_pv,   0, 1, useValues=TRUE)
    pred_npv  <- clamp(pred_npv,  0, 1, useValues=TRUE)
    pred_soil <- clamp(pred_soil, 0, 1, useValues=TRUE)
    
    # scale result by 10000
    pred_pv   <- pred_pv   * 10000
    pred_npv  <- pred_npv  * 10000
    pred_soil <- pred_soil * 10000
    
    # combine the predicted fractions
    pred_combined <- stack(pred_pv, pred_npv, pred_soil)
    # export the fraction stacks
    writeRaster(pred_combined, paste(out_dir, '/',   substring(obs_dates[obs], 2, 9),  '.tif', sep=''), datatype = "INT2S", overwrite=T)
  }
}

# =============================================================================
# list sentinel-2 rasters for unmixing
rasters_nordsea <- list.files('data/raster/sentinel-2/X0058_Y0040', full.names = TRUE)
rasters_branden <- list.files('data/raster/sentinel-2/X0068_Y0042', full.names = TRUE)
rasters_allgaeu <- list.files('data/raster/sentinel-2/X0063_Y0060', full.names = TRUE)

# =============================================================================
# apply SVRs to unmix sentinel-2 imagery, saves outpu raster for each observation date
unmix(rasters_nordsea, 'unmix/nordsea', PV_svm, NPV_svm, soil_svm)
unmix(rasters_branden, 'unmix/branden', PV_svm, NPV_svm, soil_svm)
unmix(rasters_allgaeu, 'unmix/allgaeu', PV_svm, NPV_svm, soil_svm)
