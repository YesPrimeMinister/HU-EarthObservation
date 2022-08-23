####################################################################
####################################################################

### NDFI Calculation
### Last worker: Lara
### Date: 12.07.2022

####################################################################
####################################################################


# load libraries
pck_list <- c("assertthat","GGally", "ggplot2", "tidyr", "dplyr","raster", "e1071", "rgdal",
              "viridisLite", "randomForest",  "ggpubr")
lapply(pck_list, require, character.only = TRUE)

setwd('data')
# load library endmembers
df_speclib <- read.table('spectral-library/speclib_grassland_ger.txt', header=T, sep=',')
head(df_speclib)
str(df_speclib)

# data frame changesments in order to use synthmix function
# convert class into integer (NPV as 1, PV as 2 and soil 3)
df_speclib$cover <- as.integer(as.factor(df_speclib$cover))
df_speclib

# pivot longer
# df_speclib_longer <- pivot_longer(df_speclib, cols = c(1:10), names_to = "bands")

# rename column
colnames(df_speclib)[11] <- "ClassID"
# change order of df
em_df <- df_speclib[,c(11,1:10)]
head(em_df)

###### insert shade endmember #######

synthmix <- function(df, cl_target, cl_background, n_samples=1000,
                     mix_complexity=c(2, 3, 4), p_mix_complexity=c(.7, .2, .1),
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

#### Generate synthetic mixtures
mixture_pv <- synthmix(em_df, 2, c(1,3), mix_complexity=c(2,3), n_samples = 1000,
                       p_mix_complexity=c(.5, .5), within_class_mixture=T)
mixture_npv <- synthmix(em_df, 1, c(2,3), mix_complexity=c(2,3), n_samples = 1000,
                        p_mix_complexity=c(.5, .5), within_class_mixture=T)
mixture_soil <- synthmix(em_df, 3, c(1,2), mix_complexity=c(2,3), n_samples = 1000,
                         p_mix_complexity=c(.5, .5), within_class_mixture=T)


##########################################################################################
#### Building the SVR models

# define x for training
x_pv <- subset(mixture_pv, select=-fraction)
x_npv <- subset(mixture_npv, select=-fraction)
x_soil <- subset(mixture_soil, select=-fraction)

# Suppot vector regressor training
cv <- tune.control(cross = 10)

svr_pv <- tune.svm(x=x_pv, y=mixture_pv$fraction, kernel='radial', gamma = 10^(-3:3),
                   cost = 10^(-3:3), epsilon=0.001, tunecontrol = cv)
svr_npv <- tune.svm(x=x_npv, y=mixture_npv$fraction, kernel='radial', gamma = 10^(-3:3),
                   cost = 10^(-3:3), epsilon=0.001, tunecontrol = cv)
svr_soil <- tune.svm(x=x_soil, y=mixture_soil$fraction, kernel='radial', gamma = 10^(-3:3),
                   cost = 10^(-3:3), epsilon=0.001, tunecontrol = cv)

##########################################################################################
#### Validate the models based on reference data

# load validation data
validate_path <- 'data/vector/validation/GE_VAL_s2_bands_subset.gpkg'
validate_data <- readOGR(validate_path)

# Extract validation
df_val <- data.frame(validate_data)

### extract spectral reflectances and reference fractions from dataframe
df_val_selct <- df_val[,c(2:11,30:32)]
df_val_reflectance <-
# assign names to spectral bands in validation data
names(df_val_reflectance) <- names(x_pv)
head(df_val_selct)


### apply models to the spectal reflectances
predict()
pred_svr_pv <- predict(s2, svr_pv$best.model, type='response', na.omit=T)
pred_svr_npv <- predict(s2, svr_npv$best.model, type='response', na.omit=T)
pred_svr_soil <- predict(s2, svr_soil$best.model, type='response', na.omit=T)
### compare modeled and reference fractions


##########################################################################################
#### Unmix all the available imagery

#### Unmixing the areas separately

base_dir = 'data/EO_2022_project-4_grassland-drought/raster/X0058_Y0040'
rasters_list <- list.files(paste0(base_dir, "/",pattern="*B[2-7]*.TIF$", full.names = TRUE))

unmix <- function(raster_list, models){
  
  "Function to unmix all values in a .
  
  df:                 (list) Input dataframe. First column must contain the
                      class-IDs in integer format. Remaining columns must 

  returns:            (list) Dataframe with linearily mixed features and 
                      corresponding fraction of target class (i.e. cl_target)
  "
  
  # count number of observations = number of bands in the input datasets
  n_bands <- nbands(raster_list[[1]])
  
  # init stacks for fractions at every oobservation
  stack_pv <- NaN
  stack_npv <- NaN
  stack_soil <- NaN
  
  # unmix at every time step
  for (obs in n_bands) {
    
    # read inidvidual spectral bands for individual files
    b1 <- raster(raster_list[1], band = obs)
    b2 <- raster(raster_list[2], band = obs)
    b3 <- raster(raster_list[3], band = obs)
    b4 <- raster(raster_list[4], band = obs)
    b5 <- raster(raster_list[5], band = obs)
    b6 <- raster(raster_list[6], band = obs)
    b7 <- raster(raster_list[7], band = obs)
    b8 <- raster(raster_list[8], band = obs)
    b9 <- raster(raster_list[9], band = obs)
    b10 <- raster(raster_list[10], band = obs)
    
    # stack spectral bands into a monotemporal, multispectral raster stack
    stack_monotemp <- stack(b1, b2, b3, b4, b5, b6, b7, b8, b9, b10)
    names(stack_monotemp) <- c('','')
    
    # apply unmixing models to the monotemporal raster
    pred_pv   <- predict(stack_monotemp, models[[1]], type='response', na.omit=T)
    ######## assign date as bandname - take from bandnames of original rasters
    names(pred_pv) <- ''
    pred_npv  <- predict(stack_monotemp, models[[2]], type='response', na.omit=T)
    pred_soil <- predict(stack_monotemp, models[[3]], type='response', na.omit=T)
    
    # add fractions to the overall multitemporal fraction datasets
    stack_pv <- stack(stack_pv, pred_pv)
    stack_npv <- stack(stack_npv, pred_npv)
    stack_soil <- stack(stack_soil, pred_soil)
    
  }
  return(stack_pv, stack_npv, stack_soil)
}


stack_pv <- unmix(rasters_list)


compute_ndfi <- function(stack_pv, stack_npv, stack_soil){
  
  stack_ndfi <- NaN
  
  for (obs in nbands(stack_pv)) {
    ndfi <- NaN ### add function
    names(ndfi) <- 'acquisiton name' ### add date
    stack_ndfi <- stack(stack_ndfi, ndfi)
  }
  
  
  return(stack_ndfi)
}

