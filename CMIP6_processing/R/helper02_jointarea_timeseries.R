####  Joint Survey DFO-NMFS CMIP6 Scenario Timeseries
# 9/1/2023

# Goal:
# Create mean/5th/95th timeseries information for the whole COCAFC/CROSBND
# Study Area

# This follows the workflow of helper01_regional_timeseries_construction.R






####  Packages  ####
library(here)
library(ncdf4)
library(RNetCDF)
library(janitor)
library(sf)
library(gmRi)
library(patchwork)
library(tidyverse)
library(raster)

# research path
res_path   <- cs_path("res", "")
oisst_path <- cs_path("res", "OISST/oisst_mainstays")
soda_path  <- cs_path("RES_Data", "SODA")
dfo_path <- cs_path(box_group = "mills", subfolder = "Projects/DFO_survey_data/strata_shapefiles")



####  Loading Resources  ####

####  1. Load Regions to Crop  ####

# 1. Entire VAST Area
# Just the shapefile
trawl_crsbnd <- read_sf(str_c(dfo_path, "DFO_NMFS_CRSBND_area.geojson"))





##### 1. Import/Clean Functions  ####


# load mean/5th/95th from ensembles
load_ensemble_percentiles <- function(cmip_var, ssp_scenario){
  
  # ensemble folders
  # bias corrected data folders
  ens_folders <- map(vars_neat, ~ cs_path(
    box_group = "res", 
    subfolder = str_c("CMIP6/", ssp_scenario, "/BiasCorrected/EnsembleData/", .x, "/")))
  ens_folders <- setNames(ens_folders, vars_neat)
  
  # list files:
  file_paths <- str_c(list.files(ens_folders[[cmip_var]], full.names = TRUE, pattern = "\\.grd$"))
  file_names <- list.files(ens_folders[[cmip_var]], full.names = FALSE, pattern = "\\.grd$")
  file_names <- str_remove_all(file_names, ".grd")
  file_paths <- setNames(file_paths, file_names)
  
  # load as rasters in list
  ens_stacks <- map(file_paths, raster::stack)
  
  # Return the stack
  return(ens_stacks)
}



####  Crop/Mask Functions  ####

# 1.
# Function to mask rasters using shape
mask_shape <- function(in_ras, in_mask){# Check extents
  
  # Check extent for to make sure they overlap
  # Rotate if you need to
  in_ext <- extent(in_ras)
  if(in_ext@xmax > 180){
    out_extent <- in_ext - c(360, 360, 0, 0)
    in_ras <- raster::setExtent(in_ras, out_extent)
  }
  
  # crop+mask
  r1 <- raster::crop(x = in_ras, y = in_mask)
  r2 <- raster::mask(x = r1, mask = in_mask)
  return(r2)}



# 2.
# Process means and date formats for a monthly raster stack
stack_to_df <- function(month_stack, var_name){
  var_sym <- sym(var_name)
  raster::cellStats(month_stack, "mean", na.rm = T) %>% 
    raster::as.data.frame() %>% 
    tibble::rownames_to_column() %>% 
    setNames(c("date", var_name)) %>% 
    dplyr::mutate(
      date = str_remove(date, "X"),
      date = str_replace_all(date, "[.]", "-"),
      date = as.Date(str_c(date, "-15")),
      year = lubridate::year(date))
}




####  2. Load Bias Corrected SSP Scenario Data  ####




# Load ensemble percentile data for SSP1
ssp1_sst <- load_ensemble_percentiles(cmip_var = "surf_temp", ssp_scenario = "SSP1_26")
ssp1_bt  <- load_ensemble_percentiles(cmip_var = "bot_temp", ssp_scenario = "SSP1_26")
ssp1_ss  <- load_ensemble_percentiles(cmip_var = "surf_sal", ssp_scenario = "SSP1_26")
ssp1_bs  <- load_ensemble_percentiles(cmip_var = "bot_sal", ssp_scenario = "SSP1_26")

# Load ensemble percentile data for SSP5
ssp5_sst <- load_ensemble_percentiles(cmip_var = "surf_temp", ssp_scenario = "SSP5_85")
ssp5_bt  <- load_ensemble_percentiles(cmip_var = "bot_temp", ssp_scenario = "SSP5_85")
ssp5_ss  <- load_ensemble_percentiles(cmip_var = "surf_sal", ssp_scenario = "SSP5_85")
ssp5_bs  <- load_ensemble_percentiles(cmip_var = "bot_sal", ssp_scenario = "SSP5_85")




# variables labels to pull
vars_neat <- c(
  "Surface Temperature" = "surf_temp",
  "Surface Salinity"    = "surf_sal",
  "Bottom Temperature"  = "bot_temp",
  "Bottom Salinity"     = "bot_sal")








# This was isolated to do the full combined area on its own
# Was running into namespace errors with the nesting of mask and stack steps

# 1. Mask the stacks for each variable/scenario

# SSP1
ssp1_sst_masked <- map(ssp1_sst, ~mask_shape(in_ras = .x, in_mask = trawl_crsbnd))
ssp1_bt_masked <- map(ssp1_bt, ~mask_shape(in_ras = .x, in_mask = trawl_crsbnd))
ssp1_ss_masked <- map(ssp1_ss, ~mask_shape(in_ras = .x, in_mask = trawl_crsbnd))
ssp1_bs_masked <- map(ssp1_bs, ~mask_shape(in_ras = .x, in_mask = trawl_crsbnd))

# SSP5
ssp5_sst_masked <- map(ssp5_sst, ~mask_shape(in_ras = .x, in_mask = trawl_crsbnd))
ssp5_bt_masked <- map(ssp5_bt, ~mask_shape(in_ras = .x, in_mask = trawl_crsbnd))
ssp5_ss_masked <- map(ssp5_ss, ~mask_shape(in_ras = .x, in_mask = trawl_crsbnd))
ssp5_bs_masked <- map(ssp5_bs, ~mask_shape(in_ras = .x, in_mask = trawl_crsbnd))

# # This step fails...
# crsbnd_timeseries<- purrr::map(crsbnd_masked[1], .f = stack_to_df(month_stack = .x, var_name = var_select))




# 2. Convert the masked rasters to timeseries
ssp1_sst_ts <- map_dfr(ssp1_sst_masked, ~stack_to_df(month_stack = .x, var_name = "surf_temp"), .id = "cmip_id")
ssp1_bt_ts  <- map_dfr(ssp1_bt_masked, ~stack_to_df(month_stack = .x, var_name = "bot_temp"), .id = "cmip_id")
ssp1_ss_ts  <- map_dfr(ssp1_ss_masked, ~stack_to_df(month_stack = .x, var_name = "surf_sal"), .id = "cmip_id")
ssp1_bs_ts  <- map_dfr(ssp1_bs_masked, ~stack_to_df(month_stack = .x, var_name = "bot_sal"), .id = "cmip_id")
ssp5_sst_ts <- map_dfr(ssp5_sst_masked, ~stack_to_df(month_stack = .x, var_name = "surf_temp"), .id = "cmip_id")
ssp5_bt_ts  <- map_dfr(ssp5_bt_masked, ~stack_to_df(month_stack = .x, var_name = "bot_temp"), .id = "cmip_id")
ssp5_ss_ts  <- map_dfr(ssp5_ss_masked, ~stack_to_df(month_stack = .x, var_name = "surf_sal"), .id = "cmip_id")
ssp5_bs_ts  <- map_dfr(ssp5_bs_masked, ~stack_to_df(month_stack = .x, var_name = "bot_sal"), .id = "cmip_id")



# Put all the vars in one table
crsbnd_all <- bind_rows(
  list(
    "SSP1_26" = bind_rows(
      ssp1_sst_ts,
      ssp1_bt_ts,
      ssp1_ss_ts,
      ssp1_bs_ts  
    ),
    "SSP5_85" = bind_rows(
      ssp5_sst_ts,
      ssp5_bt_ts,
      ssp5_ss_ts,
      ssp5_bs_ts
    )
    
    
  ), .id = "scenario"
)


####  Check & Save  ####

# Combine and add metadata
crsbnd_tidy <- crsbnd_all %>% 
    mutate(
      data_source = case_when(
        str_detect(cmip_id, "percentile") ~ "Ensemble Data",
        str_detect(cmip_id, "mean")       ~ "Ensemble Data",
        str_detect(cmip_id, "historic")   ~ "CMIP6 Historical",
        TRUE ~ "CMIP6 SSP"),
      ensemble_statistic = case_when(
        str_detect(cmip_id, "95th")  ~ "95th Percentile",
        str_detect(cmip_id, "5th")   ~ "5th Percentile",
        str_detect(cmip_id, "mean")  ~ "Ensemble Mean",
        TRUE ~ "Individual CMIP6 Output"),
      ensemble_statistic = factor(
        ensemble_statistic, 
        levels = c("Individual CMIP6 Output","5th Percentile",
                   "Ensemble Mean", "95th Percentile")))
  
  


# Save the file

# Folder To Put Table
jointsurvey_folder <- cs_path("res", str_c("CMIP6/", "SSP1_and_SSP5_JointSurvey_Timeseries/"))

# File Name
ts_name <- str_c(jointsurvey_folder, "CMIP6_jointsurvey_ensembles_bias_corrected.csv")
print(str_c("Saving: ", ts_name))

# Save it
write_csv(x = crsbnd_tidy, file = ts_name)



