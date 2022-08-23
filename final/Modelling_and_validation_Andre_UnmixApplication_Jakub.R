# =============================================================================
# GCG - EARTH OBSERVATION SUMMER 2022
# SESSION 08: Estimating fractional cover
# Andre Fischer
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
# load library endmembers
df_speclib <- read.table('spectral-library/speclib_grassland_ger.txt', header=T, sep=',')
head(df_speclib) 

# rescaling reflectance 
df_speclib[ ,1:10] <- df_speclib[ ,1:10] / 10000
print(df_speclib)

# add shade as a class
df_speclib[nrow(df_speclib) + 1,] <- c(0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 
                                       0.01, 0.01, 0.01, 3)
df_speclib[nrow(df_speclib) + 1,] <- c(100, 100, 100, 100, 100, 100, 100, 
                                       100, 100, 100, 3)

# Prepare data frame for synthmix-function 
df <- df_speclib
df <- df[ ,c(11, 1:10)]
colnames(df)[1] <- "class_ID"
df$class_ID <- replace(df$class_ID, df$class_ID=="PV", 1)
df$class_ID <- replace(df$class_ID, df$class_ID=="NPV", 2)
df$class_ID <- replace(df$class_ID, df$class_ID=="soil", 4)
# Now, the landcover classes are assigned to: 
# 1 = PV; 2 = NPV; 3 = shade; 4 = soil 

df$class_ID <- as.integer(df$class_ID)
# =============================================================================
# synthmix-function

synthmix <- function(df, cl_target, cl_background, n_samples,
                     mix_complexity, p_mix_complexity,
                     within_class_mixture=TRUE, include_endmember=TRUE){
  
  "Function to generate synthetic training data mixtures from pure endmember
  spectra.
  
  df:                 (list) Input dataframe. First column must contain the
                      class-IDs in integer format. Remaining columns must 
                      contain the features to be mixed. 
  cl_target:          (int) Target class' integer ID value.
  cl_background:      (int) Background class' integer ID value(s). Vector for 
                      multiple classes, e.g. 'c(2, 3, 4)'.
  n_samples:          (int) Number of synthetic training points to generate.
  mix_complexity:     (int) Vector with desired number of possible mixtures
                      between different classes.
  p_mix_complexity:   (float) Vector containing desired occurence propabilities 
                      associated to the number of possible mixtures 
                      (i.e. mix_complexity). Must be of same length as 
                      'mix_complexity' argument.
  
  returns:            (list) Dataframe with linearily mixed features and 
                      corresponding fraction of target class (i.e. cl_target)
  "
  
  # total number of classes
  all_ems <- c(cl_target, cl_background)
  n_em <- length(all_ems)
  
  # create empty df to store training data
  df_mixture <- setNames(data.frame(matrix(ncol = ncol(df), nrow = 0)), 
                         c(names(df)[2:length(df)], "fraction")) 
  
  # index list of EMs for sampling
  idx_em <- list()
  for (em in all_ems){
    idx_em[[em]] <- which(df[,1] == em)
  }
  
  # vector for fraction calculation
  zero_one <- integer(nrow(df))
  zero_one[idx_em[[cl_target]]] <- 1
  
  # iterator for generating each synthetic mixture 
  for (i in 1:n_samples) {
    
    if (length(p_mix_complexity) == 1){
      complexity = mix_complexity
    } else {
      # sample mixing complexity based on mixing likelihoods
      complexity = sample(as.vector(mix_complexity), 
                          size = 1, 
                          prob = as.vector(p_mix_complexity)) 
    }
    
    # select background EMs which will be included in the mixture
    if (within_class_mixture){
      background <- sample(all_ems, complexity - 1, replace = TRUE)
    } else {
      background <- sample(cl_background, complexity - 1, replace = FALSE)
    }
    
    # sample indices of selected EMs
    response <- c(cl_target, background)      
    drawn_index <- c()
    for (r in response){
      drawn_index <- c(drawn_index, sample(idx_em[[r]], 1))
    }
    drawn_features <- df[drawn_index, 2:length(df)]
    drawn_fraction <- zero_one[drawn_index]
    
    # sample random weights
    drawn_weights <- c()
    for (j in 1:(complexity-1)){
      if (j == 1){
        weight <- runif(1)
      } else {
        weight <- runif(1) * (1. - sum(drawn_weights))
      }
      drawn_weights <- c(drawn_weights, weight)
    }
    drawn_weights <- c(drawn_weights, (1. - sum(drawn_weights)))
    
    # calculate mixtures and associated fractions
    calc_mixtures <- apply(drawn_features * drawn_weights, 2, FUN=sum)
    calc_fraction <- sum(drawn_fraction * drawn_weights)
    
    # append to df
    df_mixture[nrow(df_mixture)+1,] <- c(calc_mixtures, calc_fraction)
  }
  
  if (include_endmember){
    df_endmember <- cbind(df[,2:length(df)], zero_one)
    colnames(df_endmember) <- c(names(df)[2:length(df)], "fraction")
    df_mixture <- rbind(df_mixture, df_endmember)
  }
  
  return(df_mixture)
  
}
# =============================================================================
# Create linear mixtures of the 3 target classes

#!! OPTION : USE 5 MIX RUNS AND CORRESPONDING MODELS AND TAKE AVERAGE

PV_mix <- synthmix(df, cl_target = 1, cl_background = c(2,3,4), n_samples = 500, 
                   mix_complexity = c(2,3,4), p_mix_complexity = c(.7, .2, .1))
NPV_mix <- synthmix(df, cl_target = 2, cl_background = c(1,3,4), n_samples = 500, 
                    mix_complexity = c(2,3,4), p_mix_complexity = c(.7, .2, .1))
soil_mix <- synthmix(df, cl_target = 4, cl_background = c(1,2,3), n_samples = 500, 
                     mix_complexity = c(2,3,4), p_mix_complexity = c(.7, .2, .1))
# =============================================================================
#### Building the SVR models
# Define accuracy from 10-fold cross-validation as optimization measure
cv <- tune.control(cross = 10) 

# tune svm and store the best model
tune_PV <- tune.svm(fraction~., data = PV_mix, kernel = 'radial', gamma = 10^(-3:3), 
                     cost = 10^(-3:3), tunecontrol = cv)
PV_svm <- tune_PV$best.model

tune_NPV <- tune.svm(fraction~., data = NPV_mix, kernel = 'radial', gamma = 10^(-3:3), 
                    cost = 10^(-3:3), tunecontrol = cv)
NPV_svm <- tune_NPV$best.model

tune_soil <- tune.svm(fraction~., data = soil_mix, kernel = 'radial', gamma = 10^(-3:3), 
                    cost = 10^(-3:3), tunecontrol = cv)
soil_svm <- tune_soil$best.model

save(PV_svm,   file='models/svr_pv_temp.txt')
save(NPV_svm,  file='models/svr_npv_temp.txt')
save(soil_svm, file='models/svr_soil_temp.txt')

load(file='models/svr_pv_500_samples.txt')
load(file='models/svr_npv_500_samples.txt')
load(file='models/svr_soil_500_samples.txt')

# =============================================================================
#### Validate the models based on reference data

# Load validation points
validation_points <- readOGR('vector/validation/GE_VAL_s2_bands_subset.gpkg')
validation_df <- as.data.frame(validation_points)
validation_df <- validation_df[ ,c(2:11,30:32)]
validation_df[1:10] <- validation_df[1:10] / 10000

# Run predictions
validation_df$PV_pred <- predict(PV_svm, validation_df[ ,1:10], na.action = na.exclude)
validation_df$NPV_pred <- predict(NPV_svm, validation_df[ ,1:10], na.action = na.exclude)
validation_df$soil_pred <- predict(soil_svm, validation_df[ ,1:10], na.action = na.exclude)

# Change values that are out of range to 0 or 1
validation_df$PV_pred[validation_df$PV_pred<0] <- 0
validation_df$PV_pred[validation_df$PV_pred>1] <- 1
validation_df$NPV_pred[validation_df$NPV_pred<0] <- 0
validation_df$NPV_pred[validation_df$NPV_pred>1] <- 1
validation_df$soil_pred[validation_df$soil_pred<0] <- 0
validation_df$soil_pred[validation_df$soil_pred>1] <- 1

# =============================================================================
# Using a linear regression to compare prediction with validation

valim_PV <- lm(validation_df$PV ~ validation_df$PV_pred, validation_df)
valim_NPV <- lm(validation_df$NPV ~ validation_df$NPV_pred, validation_df)
valim_soil <- lm(validation_df$soil ~ validation_df$soil_pred, validation_df)
# =============================================================================
# Function to plot the regressions
ggplotRegression <- function (fit, actual) {
  
  ggplot(fit$model, aes_string(x = names(fit$model)[2], y = names(fit$model)[1])) + 
    geom_point() +
    stat_smooth(method = "lm", col = "red") +
    labs(title = paste("Adj R2 = ",signif(summary(fit)$adj.r.squared, 5),
                       "Intercept =",signif(fit$coef[[1]],5 ),
                       " Slope =",signif(fit$coef[[2]], 5),
                       " P =",signif(summary(fit)$coef[2,4], 5),
                       "MAE =", mae(actual, predict(fit)), 
                       "RMSE =", rmse(actual, predict(fit))))
}
# =============================================================================
# AAAnd plotting it
ggplotRegression(valim_PV, validation_df$PV)
ggplotRegression(valim_NPV, validation_df$NPV)
ggplotRegression(valim_soil, validation_df$soil)



# =============================================================================
#### Fractional cover predictions 

unmix <- function(raster_list, bandnames, model_pv, model_npv, model_soil){
  
  "Function to unmix all values in a .
  
  df:                 (list) Input dataframe. First column must contain the
                      class-IDs in integer format. Remaining columns must 

  returns:            (list) Dataframe with linearily mixed features and 
                      corresponding fraction of target class (i.e. cl_target)
  "
  
  # count number of observations = number of bands in the input datasets
  n_bands <- nlayers(stack(raster_list[[1]]))
  obs_dates <- names(stack(raster_list[[1]]))
  
  # init stacks for fractions at every observation
  #list_pv   <- list()
  list_npv  <- list()
  #list_soil <- vector(mode='list', length=n_bands)
  
  print(paste('Processing observation # out of ', n_bands, ':', sep=''))
  # unmix at every time step
  for (obs in 1:n_bands) {
    
    print(paste(obs, '/', n_bands, sep=''))
    # read inidvidual spectral bands for individual files
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
    names(stack_monotemp) <- bandnames
    
    # apply unmixing models to the monotemporal raster
    #pred_pv   <- predict(stack_monotemp, model_pv, type='response', na.omit=T)
    pred_npv  <- predict(stack_monotemp, model_npv, type='response', na.omit=T)
    #pred_soil <- predict(stack_monotemp, model_soil, type='response', na.omit=T)
    
    # add raster band name according to date
    #names(pred_pv)   <- substring(obs_dates[obs], 1, 9)
    names(pred_npv)  <- substring(obs_dates[obs], 1, 9)
    #names(pred_soil) <- substring(obs_dates[obs], 1, 9)
    
    print(pred_npv)
    # remove values above 1 / below 0
    #pred_pv   <- clamp(pred_pv,   0, 1, useValues=TRUE)
    pred_npv  <- clamp(pred_npv,  0, 1, useValues=TRUE)
    #pred_soil <- clamp(pred_soil, 0, 1, useValues=TRUE)
    
    # scale result by 10000
    #pred_pv   <- pred_pv   * 10000
    pred_npv  <- pred_npv  * 10000
    #pred_soil <- pred_soil * 10000
    
    # add fractions to the overall multitemporal fraction datasets
    #list_pv[obs] <- pred_pv
    list_npv[obs] <- pred_npv
    #list_soil[obs] <- pred_soil
  }
  
  # stack all the individually created layers
  #stack_pv   <- stack(list_pv)
  stack_npv  <- stack(list_npv)
  #stack_soil <- stack(list_soil)
  
  # convert stacks to terra rast for correct export - otherwise band names are not preserved
  #stack_pv   <- rast(list_pv)
  #stack_npv  <- sds(stack_npv)
  #stack_soil <- sds(stack_soil)
  
  return(stack_pv)
  #unmixed_stacks <- c(stack_pv, stack_npv, stack_soil)
  
  #return(unmixed_stacks)
}



# HERE; WE CAN LOAD THE SENTINEL 2 IMAGERY, THAT WE WILL DECIDE TO USE
TEMP_dir = 'E:/__Earth_Observation/grasslands/data/raster/sentinel-2/X0058_Y0040'
rasters_list <- list.files(TEMP_dir, full.names = TRUE)
rasters_list
band_names <- dimnames(PV_svm$SV)[[2]]

# testing unmix
TEMP_unmixed <- unmix(rasters_list, band_names, PV_svm, NPV_svm, soil_svm)

writeRaster(TEMP_unmixed, "E:/__Earth_Observation/grasslands/data/nordsea_npv_500.tif", datatype = "INT2S")

TEMP_unmixed_terra <- unmix_terra(rasters_list, band_names, PV_svm, NPV_svm, soil_svm)
TEMP_unmixed_individual <- unmix_individual(rasters_list, band_names, model_npv=NPV_svm)




TEMP_terra <- rast("E:/__Earth_Observation/grasslands/data/nordsea_npv_500.tif")
obs_names <- rast(rasters_list[1])
names(TEMP_terra) <- obs_names@ptr$names
writeRaster(TEMP_terra, "E:/__Earth_Observation/grasslands/data/nordsea_npv_500_terra.tif", datatype = "INT2S")

# extract spectral band names from a svr model
band_names <- dimnames(PV_svm$SV)[[2]]

# List sentinel-2 filenames for individual areas
rasters_nordsea <- list.files('data/raster/sentinel-2/X0058_Y0040', full.names = TRUE)
rasters_branden <- list.files('data/raster/sentinel-2/X0068_Y0042', full.names = TRUE)
rasters_allgaeu <- list.files('data/raster/sentinel-2/X0063_Y0060', full.names = TRUE)

# apply SVRs to unmix sentinel-2 imagery
unmixed_nordsea <- unmix(rasters_nordsea, band_names, PV_svm, NPV_svm, soil_svm)
unmixed_branden <- unmix(rasters_branden, band_names, PV_svm, NPV_svm, soil_svm)
unmixed_allgaeu <- unmix(rasters_allgaeu, band_names, PV_svm, NPV_svm, soil_svm)

writeRaster(tcb_sd_raster, "tcb_sd_08_11.tif", datatype = "INT2S",overwrite = T )


unmix_individual <- function(raster_list, bandnames, model_pv=F, model_npv=F, model_soil=F){
  
  "Function to unmix all values in a .
  
  df:                 (list) Input dataframe. First column must contain the
                      class-IDs in integer format. Remaining columns must 

  returns:            (list) Dataframe with linearily mixed features and 
                      corresponding fraction of target class (i.e. cl_target)
  "
  
  # count number of observations = number of bands in the input datasets
  n_bands <- nlayers(stack(raster_list[[1]]))
  obs_dates <- names(stack(raster_list[[1]]))
  
  # init stacks for fractions at every observation
  if (model_pv   != F) {list_pv   <- list()}
  if (model_npv  != F) {list_npv  <- list()}
  if (model_soil != F) {list_soil <- list()}
  
  print(paste('Processing observation # out of ', n_bands, ':', sep=''))
  
  # unmix at every time step
  for (obs in 1:n_bands) {
    
    print(obs)
    # read individual spectral bands for individual files
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
    names(stack_monotemp) <- bandnames
    
    
    if (model_pv != F) {
      pred_pv   <- predict(stack_monotemp, model_pv, type='response', na.omit=T)
      names(pred_pv)   <- substring(obs_dates[obs], 1, 9)
      
      pred_pv   <- clamp(pred_pv,   0, 1, useValues=TRUE)
      pred_pv   <- pred_pv   * 10000
      
      list_pv[obs] <- pred_pv
    }
    
    if (model_npv != F) {
      pred_npv  <- predict(stack_monotemp, model_npv, type='response', na.omit=T)
      names(pred_npv)  <- substring(obs_dates[obs], 1, 9)
      
      pred_npv  <- clamp(pred_npv,  0, 1, useValues=TRUE)
      pred_npv  <- pred_npv  * 10000
      
      list_npv[obs] <- pred_npv
    }
    
    if (model_soil != F) {
      pred_soil <- predict(stack_monotemp, model_soil, type='response', na.omit=T)
      names(pred_soil) <- substring(obs_dates[obs], 1, 9)
      
      pred_soil <- clamp(pred_soil, 0, 1, useValues=TRUE)
      pred_soil <- pred_soil * 10000
      
      list_soil[obs] <- pred_soil
    }
    # apply unmixing models to the monotemporal raster
    # add raster band name according to date
    # remove values above 1 / below 0
    # scale result by 10000
    # add fractions to the overall multitemporal fraction datasets
  }
  
  
  # stack all the individually created layers
  if (model_pv   != F) {stack_pv   <- stack(list_pv)}
  if (model_npv  != F) {stack_npv  <- stack(list_npv)}
  if (model_soil != F) {stack_soil <- stack(list_soil)}
  
  # convert stacks to terra rast for correct export - otherwise band names are not preserved
  #if (model_pv)   {stack_pv   <- sds(stack_pv)}
  #if (model_npv)  {stack_npv  <- sds(stack_npv)}
  #if (model_soil) {stack_soil <- sds(stack_soil)}
  
  
  unmixed_stacks <- c()
  if (model_pv   != F) {append(unmixed_stacks, stack_pv)}
  if (model_npv  != F) {append(unmixed_stacks, stack_npv)}
  if (model_soil != F) {append(unmixed_stacks, stack_soil)}
  
  return(unmixed_stacks)
}



unmix_terra <- function(raster_list, bandnames, model_pv, model_npv, model_soil){
  
  "Function to unmix all values in a .
  
  df:                 (list) Input dataframe. First column must contain the
                      class-IDs in integer format. Remaining columns must 

  returns:            (list) Dataframe with linearily mixed features and 
                      corresponding fraction of target class (i.e. cl_target)
  "
  
  # count number of observations = number of bands in the input datasets
  n_bands <- nlayers(stack(raster_list[[1]]))
  obs_dates <- names(stack(raster_list[[1]]))
  
  # init stacks for fractions at every observation
  list_pv   <- vector(mode='list', length=n_bands)
  #list_npv  <- vector(mode='list', length=n_bands)
  #list_soil <- vector(mode='list', length=n_bands)
  
  print(paste('Processing observation # out of ', n_bands, ':', sep=''))
  # unmix at every time step
  for (obs in 1:n_bands) {
    
    print(paste(obs, '/', n_bands, sep=''))
    # read inidvidual spectral bands for individual files
    blue  <- rast(raster_list[1], lyrs = obs)
    green <- rast(raster_list[3], lyrs = obs)
    red   <- rast(raster_list[8], lyrs = obs)
    re1   <- rast(raster_list[5], lyrs = obs)
    re2   <- rast(raster_list[6], lyrs = obs)
    re3   <- rast(raster_list[7], lyrs = obs)
    bnir  <- rast(raster_list[2], lyrs = obs)
    nir   <- rast(raster_list[4], lyrs = obs)
    swir1 <- rast(raster_list[9], lyrs = obs)
    swir2 <- rast(raster_list[10], lyrs = obs)
    
    # stack spectral bands into a monotemporal, multispectral raster stack
    stack_monotemp <- c(blue, green, red, re1, re2, re3, bnir, nir, swir1, swir2)
    names(stack_monotemp) <- bandnames
    
    # apply unmixing models to the monotemporal raster
    pred_pv   <- predict(stack_monotemp, model_pv, type='response', na.omit=T)
    #pred_npv  <- predict(stack_monotemp, model_npv, type='response', na.omit=T)
    #pred_soil <- predict(stack_monotemp, model_soil, type='response', na.omit=T)
    
    # add raster band name according to date
    names(pred_pv)   <- substring(obs_dates[obs], 1, 9)
    #names(pred_npv)  <- substring(obs_dates[obs], 1, 9)
    #names(pred_soil) <- substring(obs_dates[obs], 1, 9)
    
    # remove values above 1 / below 0
    pred_pv   <- clamp(pred_pv,   0, 1, values=TRUE)
    #pred_npv  <- clamp(pred_npv,  0, 1, useValues=TRUE)
    #pred_soil <- clamp(pred_soil, 0, 1, useValues=TRUE)
    
    # scale result by 10000
    pred_pv   <- pred_pv   * 10000
    #pred_npv  <- pred_npv  * 10000
    #pred_soil <- pred_soil * 10000
    
    #pred_pv <- rast(pred_pv)
    
    # add fractions to the overall multitemporal fraction datasets
    list_pv[obs] <- pred_pv
    #append(list_pv, pred_pv)
    #list_npv[obs] <- pred_npv
    #list_soil[obs] <- pred_soil
  }
  
  # stack all the individually created layers
  #stack_pv   <- stack(list_pv)
  #stack_npv  <- stack(list_npv)
  #stack_soil <- stack(list_soil)
  
  # convert stacks to terra rast for correct export - otherwise band names are not preserved
  stack_pv   <- rast(list_pv)
  #stack_npv  <- sds(stack_npv)
  #stack_soil <- sds(stack_soil)
  
  return(stack_pv)
  #unmixed_stacks <- c(stack_pv, stack_npv, stack_soil)
  
  #return(unmixed_stacks)
}



raster_list <- rasters_list
obs = 27
model_pv <- PV_svm
model_npv <- NPV_svm
model_soil <- soil_svm

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
pred_pv   <- predict(stack_monotemp, model_pv, type='response', na.omit=T)
pred_npv  <- predict(stack_monotemp, model_npv, type='response', na.omit=T)
pred_soil <- predict(stack_monotemp, model_soil, type='response', na.omit=T)

# add raster band name according to date
#names(pred_pv)   <- substring(obs_dates[obs], 1, 9)
#names(pred_npv)  <- substring(obs_dates[obs], 1, 9)
#names(pred_soil) <- substring(obs_dates[obs], 1, 9)

print(pred_npv)
# remove values above 1 / below 0
pred_pv   <- clamp(pred_pv,   0, 1, useValues=TRUE)
pred_npv  <- clamp(pred_npv,  0, 1, useValues=TRUE)
pred_soil <- clamp(pred_soil, 0, 1, useValues=TRUE)

# scale result by 10000
pred_pv   <- pred_pv   * 10000
pred_npv  <- pred_npv  * 10000
pred_soil <- pred_soil * 10000

tempstack <- stack(pred_pv, pred_npv, pred_soil)

writeRaster(tempstack, "E:/__Earth_Observation/grasslands/data/nordsea_20180724_500samples.tif", datatype = "INT2S")



