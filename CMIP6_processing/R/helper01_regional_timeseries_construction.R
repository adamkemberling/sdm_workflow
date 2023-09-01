####  ABOUT:  ####
####  Bias-Corrected Spaghetti Plots
####  11/5/2021



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



####  Loading Resources  ####

####  1. Load Regions to Crop  ####


# 1. Gulf of Maine
trawl_regions <- get_timeseries_paths("nmfs_trawl_regions", box_location = "cloudstorage")
trawl_gom <- read_sf(trawl_regions$gulf_of_maine$shape_path)


# 2. The Full survey area we use
# Load all the strata and just filter out the crap ones
trawl_full <- read_sf(str_c(res_path, "Shapefiles/BottomTrawlStrata/BTS_Strata.shp"))  %>% 
  janitor::clean_names() %>% 
  filter(strata >= 01010 ,
         strata <= 01760,
         strata != 1310,
         strata != 1320,
         strata != 1330,
         strata != 1350,
         strata != 1410,
         strata != 1420,
         strata != 1490) 


# 3. DFO Survey Area
dfo_path <- cs_path(box_group = "mills", subfolder = "Projects/DFO_survey_data/strata_shapefiles")
dfo_area <- read_sf(str_c(dfo_path, "MaritimesRegionEcosystemAssessmentBoundary.shp"))


# 4. Ecological Production Units
epu_sf <- ecodata::epu_sf



# 5. Entire VAST Area

# Just the shapefile
trawl_crsbnd <- read_sf(str_c(dfo_path, "DFO_NMFS_CRSBND_area.geojson"))








##### Import Functions  ####

# loading specific variables
load_bias_corrected <- function(cmip_var, ssp_scenario){
  
  # bias corrected data folders
  var_folders <- map(vars_neat, ~ cs_path(
    box_group = "res", 
    subfolder = str_c("CMIP6/", ssp_scenario, "/BiasCorrected/IndividualModels/", .x, "/")))
  var_folders <- setNames(var_folders, vars_neat)
  
  
  # list files:
  file_paths <- str_c(list.files(var_folders[[cmip_var]], full.names = TRUE, pattern = "\\.grd$"))
  file_names <- list.files(var_folders[[cmip_var]], full.names = FALSE, pattern = "\\.grd$")
  file_names <- str_remove_all(file_names, ".grd")
  file_paths <- setNames(file_paths, file_names)
  
  # load as rasters in list
  var_stacks <- map(file_paths, raster::stack)
  
  # Return the stack
  return(var_stacks)
}


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



####  3. Set Options  ####

# This code is run sequentially for: SSP scenario & ssp var



####  Select SSP Scenario  ####
# ssp_select <- "SSP1_26"
ssp_select <- "SSP5_85"



####  Select Variable  ####

# Select ONE variable to use for workflow:
# var_select <- "surf_temp"
# var_select <- "surf_sal"
# var_select <- "bot_temp"
var_select <- "bot_sal"



# variables to pull
vars_neat <- c(
  "Surface Temperature" = "surf_temp",
  "Surface Salinity"    = "surf_sal", 
  "Bottom Temperature"  = "bot_temp",
  "Bottom Salinity"     = "bot_sal")




####  Load Bias Corrected SSP Scenario Data  ####


# load bias corrected data for individual scenario runs
observed_var <- load_bias_corrected(cmip_var = var_select, ssp_scenario = ssp_select)

# Load ensemble percentile data
ensemble_var <- load_ensemble_percentiles(cmip_var = var_select, ssp_scenario = ssp_select)






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





####### Singe Region Timeseries Debugging  ##########

# This was isolated to do the full combined area on its own
# Was running into namespace errors with the nesting of mask and stack steps

# Run a single region and skip the preamble problems
crsbnd_scenarios <- map(observed_var, ~mask_shape(in_ras = .x, in_mask = trawl_crsbnd))
crsbnd_percentiles <- map(observed_var, ~mask_shape(in_ras = .x, in_mask = trawl_crsbnd))

# # This step fails...
# crsbnd_timeseries<- purrr::map(crsbnd_masked[1], .f = stack_to_df(month_stack = .x, var_name = var_select))

# Works exactly fine outside of the function...
scenarios_test <- map_dfr(crsbnd_masked, ~raster::cellStats(.x, "mean", na.rm = T) %>% 
                 as.data.frame()%>% 
                 tibble::rownames_to_column() %>% 
                 setNames(c("date", var_select)) %>% 
                 dplyr::mutate(
                   date = str_remove(date, "X"),
                   date = str_replace_all(date, "[.]", "-"),
                   date = as.Date(str_c(date, "-15")),
                   year = lubridate::year(date)),
                 .id = "cmip_id") %>% 
  mutate(Region = "Joint NMFS-DFO Survey Area")






################################
################################

#################################
####  Perform the Data Prep for All Regions  ####


# 3. Function to combine steps and run for all areas
# Mask and Stack together
# put them together, set which regions are going down the pipeline here
mask_and_stack <- function(masking_var, var_name){
  
  # Mask and get timeseries for each area
  
  # Gulf of Maine
  masked_gom   <- mask_shape(in_ras = masking_var, in_mask =  trawl_gom)
  masked_gom   <- stack_to_df(month_stack = masked_gom, var_name =  var_name)
  
  # Trawl Survey
  masked_trawl <- mask_shape(in_ras = masking_var, in_mask =  trawl_full)
  masked_trawl <- stack_to_df(month_stack = masked_trawl,var_name =  var_name)
  
  # DFO Survey
  masked_dfo   <- mask_shape(in_ras = masking_var, in_mask =  dfo_area)
  masked_dfo   <- stack_to_df(month_stack = masked_dfo, var_name = var_name)
  
  # EPUs
  masked_gom_epu <- mask_shape(in_ras = masking_var, epu_sf %>% filter(EPU == "GOM"))
  masked_gom_epu <- stack_to_df(month_stack = masked_gom_epu, var_name = var_name)
  masked_gb_epu  <- mask_shape(in_ras = masking_var, epu_sf %>% filter(EPU == "GB"))
  masked_gb_epu  <- stack_to_df(month_stack = masked_gb_epu, var_name = var_name)
  masked_ss_epu  <- mask_shape(in_ras = masking_var, epu_sf %>% filter(EPU == "SS"))
  masked_ss_epu  <- stack_to_df(month_stack = masked_ss_epu, var_name = var_name)
  masked_mab_epu <- mask_shape(in_ras = masking_var, epu_sf %>% filter(EPU == "MAB"))
  masked_mab_epu <- stack_to_df(month_stack = masked_mab_epu, var_name = var_name)
  
  
  # FULL DFO + NMFS Region
  masked_crsbnd <- mask_shape(masking_var, trawl_crsbnd)
  masked_crsbnd <- stack_to_df(masked_mab_epu, var_name)
  
  
  # Put in list, combine
  bind_rows(list(
    "Gulf of Maine"          = masked_gom,
    "US Survey Area"         = masked_trawl,
    "Canadian Survey Area"   = masked_dfo,
    "EPU_GOM"                = masked_gom_epu,
    "EPU_GB"                 = masked_gb_epu,
    "EPU_SS"                 = masked_ss_epu,
    "EPU_MAB"                = masked_mab_epu,
    "combined_surveys"       = masked_crsbnd
  ), .id = "Region")
}





# Runs one variable at a time:
masked_var <- map_dfr(
  .x = observed_var,
  .f = ~ mask_and_stack(masking_var = .x, var_name = var_select), 
  .id = "cmip_id")


# Run the percentiles too
masked_percentiles <- map_dfr(
  .x = ensemble_var,
  .f = ~ mask_and_stack(masking_var = .x, var_name = var_select), 
  .id = "cmip_id")





####  Check & Save  ####

# Combine and add metadata
var_combined <- bind_rows(
  list(masked_var, masked_percentiles)) %>% 
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
                 "Ensemble Mean", "95th Percentile"))
    )



##### Check With Plot  ####
# palette options
spaghetti_pal <- c(
  "Individual CMIP6 Output" = "gray",
  "5th Percentile" = "lightblue",
  "Ensemble Mean" = as.character(gmri_cols("gmri blue")),
  "95th Percentile" = "dark red")

spaghetti_sizes <- c(
  "Individual CMIP6 Output" = 0.3,
  "5th Percentile" = .75,
  "Ensemble Mean" = .75,
  "95th Percentile" = .75)

spaghetti_alpha <- c(
  "Individual CMIP6 Output" = 0.2,
  "5th Percentile" = .6,
  "Ensemble Mean" = .6,
  "95th Percentile" = .6)

# format variable for display on axis
var_titles <- c(
  "bot_temp" = expression("Bottom Temperature"~degree~"C"),
  "surf_temp" = expression("Surface Temperature"~degree~"C"),
  "bot_sal" = "Bottom Salinity",
  "surf_sal" = "Surface Salinity")
var_label <- var_titles[var_select]


# And Plot!
ggplot(data = filter(var_combined, ensemble_statistic == "Individual CMIP6 Output"), 
       aes_string("date", var_select, group = "cmip_id", 
                  color = "ensemble_statistic",
                  size = "ensemble_statistic",
                  alpha = "ensemble_statistic")) +
  geom_line() +
  geom_line(data = filter(var_combined, ensemble_statistic != "Individual CMIP6 Output")) +
  scale_color_manual(values = spaghetti_pal) +
  scale_size_manual(values = spaghetti_sizes) +
  scale_alpha_manual(values = spaghetti_alpha) +
  facet_wrap(~Region, ncol = 1, scales = "free") +
  labs(x = "", y = var_label, color = "Ensemble Statistic") +
  guides(
    color = guide_legend(title.position = "top", title.hjust = 0.5, nrow = 2),
    alpha = "none",
    size = "none")





####  Export for later:  ####

# Folder To Put Table
scenario_folder <- cs_path("res", str_c("CMIP6/", ssp_select, "/BiasCorrected/TimeseriesData/"))

# File Name
ts_name <- str_c(scenario_folder, "CMIP6_bias_corrected_regional_", var_select, ".csv")
print(str_c("Saving: ", ts_name))

# Save it
write_csv(x = var_combined, file = ts_name)






####______________________####


# Prepare the reference datasets the same way so they can be compared against them


####  Reference Datasets  ####

# Mask the obersvarional datasets as their own timeseries 
# so we can compare quickly without masking etc.


# SODA vars
soda_ssal  <- stack(str_c(soda_path, "SODA_Salt_Red.nc")) # Has 50 depths, defaults to level 1 (surface)
soda_bsal  <- stack(str_c(soda_path, "SODA_Salt_Red_bottomLayer.nc"))
soda_btemp <- stack(str_c(soda_path, "SODA_Temp_Red_bottomLayer.nc"))


# Load Monthly data
oisst_month_path <- cs_path("res", "OISST/oisst_mainstays/monthly_averages")
oisst_monthly <- stack(str_c(oisst_month_path, "oisst_monthly.nc"), varname = "sst")





# Perform the masking:

# Runs one variable at a time:
masked_surface_sal <- mask_and_stack(masking_var = soda_ssal, var_name = "surf_sal")
masked_bottom_sal  <- mask_and_stack(masking_var = soda_bsal, var_name = "bot_sal")
masked_bottom_temp <- mask_and_stack(masking_var = soda_btemp, var_name = "bot_temp")
masked_surf_temp   <- mask_and_stack(masking_var = oisst_monthly, var_name = "surf_temp")


# make months standardize
tune_months <- function(x){
  x  %>% 
    mutate(month =  str_pad(month(date), side = "left", pad = "0", width = 2), 
           date = as.Date(str_c(year,"-", month, "-15")))
}

masked_ssal <- tune_months(masked_surface_sal)
masked_stemp <- tune_months(masked_surf_temp)
masked_bsal <- tune_months(masked_bottom_sal)
masked_btemp <- tune_months(masked_bottom_temp)


####  Check & Save  ####

# Combine and add metadata
# "date" is going to need some tuning before join
references_combined <- left_join(
  masked_ssal, masked_stemp, join_by("Region", "date", "year", "month")) %>% 
  left_join(masked_bsal,join_by("Region", "date", "year", "month")) %>% 
  left_join(masked_btemp, join_by("Region", "date", "year", "month"))



# Save these things out
references_folder <- str_c(res_path, "CMIP6/Bias_Correction_Checking/bias_correction_reference_dataset_monthly_regional_means.csv")
write_csv(
  x = references_combined, 
  file = references_folder
)

