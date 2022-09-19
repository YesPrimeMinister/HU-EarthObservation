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
# set directories

# Set working directory
dir_main <- 'data'
setwd(dir_main)

# input raster dirs
unmixed_nordsea <- 'data/raster/unmix/nordsea'
unmixed_branden <- 'data/raster/unmix/branden'
unmixed_allgaeu <- 'data/raster/unmix/allgaeu'

# output ndfi dirs
ndfi_nordsea <- 'data/raster/ndfi/nordsea'
ndfi_branden <- 'data/raster/ndfi/branden'
ndfi_allgaeu <- 'data/raster/ndfi/allgaeu'

# output ndfi stack path
ndfi_nordsea_stack <- 'data/raster/ndfi/ndfi_nordsea.tif'
ndfi_branden_stack <- 'data/raster/ndfi/ndfi_branden.tif'
ndfi_allgaeu_stack <- 'data/raster/ndfi/ndfi_allgaeu.tif'


unmixed_nordsea <- 'E:/__Earth_Observation/grasslands/data/raster/unmix/nordsea'
unmixed_branden <- 'E:/__Earth_Observation/grasslands/data/raster/unmix/branden'
unmixed_allgaeu <- 'E:/__Earth_Observation/grasslands/data/raster/unmix/allgaeu'
ndfi_nordsea <- 'E:/__Earth_Observation/grasslands/data/raster/ndfi/nordsea'
ndfi_branden <- 'E:/__Earth_Observation/grasslands/data/raster/ndfi/branden'
ndfi_allgaeu <- 'E:/__Earth_Observation/grasslands/data/raster/ndfi/allgaeu'
ndfi_nordsea_stack <- 'E:/__Earth_Observation/grasslands/data/raster/ndfi/ndfi_nordsea.tif'
ndfi_branden_stack <- 'E:/__Earth_Observation/grasslands/data/raster/ndfi/ndfi_branden.tif'
ndfi_allgaeu_stack <- 'E:/__Earth_Observation/grasslands/data/raster/ndfi/ndfi_allgaeu.tif'


# =============================================================================
# compute minimum in a period between April and half of July

compute_minimum_npv <- function(in_dir) {
  # load file names for individual months
  paths_april     <- list.files(in_dir, pattern='201804', full.names=T)
  paths_may       <- list.files(in_dir, pattern='201805', full.names=T)
  paths_june_half <- list.files(in_dir, pattern='201806[01:15]', full.names=T)
  # combine file names into one vector
  paths <- c(paths_april, paths_may, paths_june_half)
  
  # create empty list to hold loaded rasters
  raster_list <- list()
  # load individual rasters
  for (obs in 1:length(paths)) {
    loaded_raster <- raster(paths[obs], band = 2)
    raster_list[obs] <- loaded_raster
  }
  # stack the loaded rasters
  raster_stack <- stack(raster_list)
  # get minimum value from raster stack
  raster_min <- min(raster_stack, na.rm=T)
  
  return(raster_min)
}

min_npv_nordsea <- compute_minimum_npv(unmixed_nordsea)
min_npv_branden <- compute_minimum_npv(unmixed_branden)
min_npv_allgaeu <- compute_minimum_npv(unmixed_allgaeu)

# =============================================================================
# compute ndfi

compute_ndfi <- function(in_dir, out_dir, min_npv) {
  in_paths <- list.files(in_dir, full.names=T)
  total_obs <- length(in_paths)

  print(paste('Processing observation # out of ', total_obs, ':', sep=''))
  
  for (idx in 1:total_obs) {
    print(paste(idx, '/', total_obs, sep=''))
    
    # load raster and observation date
    x <- stack(in_paths[idx])
    obs <- substring(names(x)[1], 2, 9)
    
    # scale in raster and rename layers
    min_npv_flt <- min_npv / 10000
    x           <- x       / 10000
    names(x) <- c('pv', 'npv', 'soil')
    
    # compute ndfi
    ndfi <- (((x$npv - min_npv_flt) + x$soil) - x$pv) / ((x$npv - min_npv_flt) + x$pv + x$soil)
    
    # scale and export the resulting ndfi raster
    ndfi <- ndfi * 10000
    writeRaster(ndfi, paste(out_dir, '/', obs,  '.tif', sep=''), datatype = "INT2S", overwrite=T)
    
  }
}

compute_ndfi(unmixed_nordsea, ndfi_nordsea, min_npv_nordsea)
compute_ndfi(unmixed_branden, ndfi_branden, min_npv_branden)
compute_ndfi(unmixed_allgaeu, ndfi_allgaeu, min_npv_allgaeu)

# =============================================================================
# stack resulting rasters

stack_ndfi <- function(in_dir, out_path) {
  
  paths  <- list.files(in_dir, full.names=T)
  fnames <- list.files(in_dir, full.names=F)
  raster_list <- list()
  
  # load individual rasters
  for (obs in 1:length(paths)) {
    loaded_raster <- rast(paths[obs])
    names(loaded_raster) <- paste('ndfi_', substring(fnames[obs], 1, 8), sep='')
    raster_list[obs] <- loaded_raster
  }
  # stack the loaded rasters
  raster_stack <- rast(raster_list)
  writeRaster(raster_stack, out_path, datatype = "INT2S", overwrite=T)
}

stack_ndfi(ndfi_nordsea, ndfi_nordsea_stack)
stack_ndfi(ndfi_branden, ndfi_branden_stack)
stack_ndfi(ndfi_allgaeu, ndfi_allgaeu_stack)

