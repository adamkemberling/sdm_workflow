####  ABOUT:  ####
####  Bias-Corrected Spaghetti Plots
####  11/5/2021



####  Packages  ####
library(here)
library(ncdf4)
library(RNetCDF)
library(raster)
library(janitor)
library(sf)
library(gmRi)
library(patchwork)
library(tidyverse)

# research path
res_path <- box_path("res", "")

# Set theme for plots 
theme_set(theme_bw() +
            theme(legend.position = "bottom",
                  strip.text = element_text(color = "white", 
                                            face = "bold"),
                  strip.background = element_rect(
                    color = "white", 
                    fill = "#36454F", 
                    size = 0.75, 
                    linetype="solid"),
                  legend.background = element_rect(fill = "transparent", 
                                                   color = "black")))



####  Select Variable  ####

# Select ONE variable to use for workflow:
# var_select <- "surf_temp"
# var_select <- "surf_sal"
# var_select <- "bot_temp"
var_select <- "bot_sal"





####  Load Bias Corrected Data  ####

# variables to pull
vars_neat <- c(
   "Surface Temperature" = "surf_temp",
   "Surface Salinity" = "surf_sal", 
   "Bottom Temperature" = "bot_temp",
   "Bottom Salinity" = "bot_sal")



##### Import Functions  ####

# load variable
load_bias_corrected <- function(cmip_var){
  
  # bias corrected data folders
  var_folders <- map(vars_neat, ~ box_path("res", str_c("CMIP6/SSP5_85/BiasCorrected/IndividualModels/", .x, "/")))
  var_folders <- setNames(var_folders, vars_neat)
  
  
  # list files:
  file_paths <- str_c(list.files(var_folders[[cmip_var]], full.names = TRUE, pattern = "\\.grd$"))
  file_names <- list.files(var_folders[[cmip_var]], full.names = FALSE, pattern = "\\.grd$")
  file_names <- str_remove_all(file_names, ".grd")
  file_paths <- setNames(file_paths, file_names)
  
  # load as rasters in list
  var_stacks <- map(file_paths, raster::stack)
  
  return(var_stacks)
}


# load ensemble data
load_ensemble_percentiles <- function(cmip_var){
  
  # ensemble folders
  # bias corrected data folders
  ens_folders <- map(vars_neat, ~ box_path("res", str_c("CMIP6/SSP5_85/BiasCorrected/EnsembleData/", .x, "/")))
  ens_folders <- setNames(ens_folders, vars_neat)
  
  # list files:
  file_paths <- str_c(list.files(ens_folders[[cmip_var]], full.names = TRUE, pattern = "\\.grd$"))
  file_names <- list.files(ens_folders[[cmip_var]], full.names = FALSE, pattern = "\\.grd$")
  file_names <- str_remove_all(file_names, ".grd")
  file_paths <- setNames(file_paths, file_names)
  
  # load as rasters in list
  ens_stacks <- map(file_paths, raster::stack)
  
  return(ens_stacks)
}




##### Import data for variable

# load bias corrected data
observed_var <- load_bias_corrected(cmip_var = var_select)

# Load ensemble percentile data
ensemble_var <- load_ensemble_percentiles(cmip_var = var_select)


####  Crop for Each Region  ####

# Getting Timeseries and Shapefile Locations via gmRi
# returns path to oisst timeseries & path to shapefile:
trawl_regions <- get_timeseries_paths("nmfs_trawl_regions")

# # NMFS Trawl Regions
trawl_gom <- read_sf(trawl_regions$gulf_of_maine$shape_path)


# Load all the strata and just filter out the crap ones
trawl_full <- read_sf(str_c(res_path, "Shapefiles/BottomTrawlStrata/BTS_Strata.shp"))  %>% 
  clean_names() %>% 
  filter(strata >= 01010 ,
         strata <= 01760,
         strata != 1310,
         strata != 1320,
         strata != 1330,
         strata != 1350,
         strata != 1410,
         strata != 1420,
         strata != 1490) 


# DFO Data
dfo_path <- shared.path(group = "Mills Lab", folder = "Projects/DFO_survey_data/strata_shapefiles")
dfo_area <- read_sf(str_c(dfo_path, "MaritimesRegionEcosystemAssessmentBoundary.shp"))



# Function to mask rasters using shape
mask_shape <- function(in_ras, in_mask){# Check extents
  # Check extent for to make sure they overlap
  in_ext <- extent(in_ras)
  if(in_ext@xmax > 180){
    out_extent <- in_ext - c(360, 360, 0, 0)
    in_ras <- setExtent(in_ras, out_extent)
  }
  
  # crop+mask
  r1 <- raster::crop(x = in_ras, y = in_mask)
  r2 <- raster::mask(x = r1, mask = in_mask)
  return(r2)}


# Process means and date formats for a monthly raster stack
stack_to_df <- function(month_stack, var_name){
  var_sym <- sym(var_name)
  cellStats(month_stack, mean, na.rm = T) %>% 
    as.data.frame() %>% 
    rownames_to_column() %>% 
    setNames(c("date", var_name)) %>% 
    mutate(date = str_remove(date, "X"),
           date = str_replace_all(date, "[.]", "-"),
           date = as.Date(str_c(date, "-15")),
           year = lubridate::year(date)) 
}







####  Mask and Get Timeseries  ####

# Runs one variable at a time:
masked_var <- map_dfr(observed_var, function(masking_var, var_name = var_select){
  
  # Mask and get timeseries for each area
  masked_gom   <- mask_shape(masking_var, trawl_gom)
  masked_gom   <- stack_to_df(masked_gom, var_name)
  masked_trawl <- mask_shape(masking_var, trawl_full)
  masked_trawl <- stack_to_df(masked_trawl, var_name)
  masked_dfo   <- mask_shape(masking_var, dfo_area)
  masked_dfo   <- stack_to_df(masked_dfo, var_name)
  
  # Put in list, combine
  bind_rows(list(
    "Gulf of Maine"  = masked_gom,
    "US Survey Area" = masked_trawl,
    "Canadian Survey Area"   = masked_dfo), .id = "Region")
  
}, .id = "cmip_id")


# Run the percentiles too
masked_percentiles <- map_dfr(ensemble_var, function(masking_var, var_name = var_select){
  
  # Mask and get timeseries for each area
  masked_gom   <- mask_shape(masking_var, trawl_gom)
  masked_gom   <- stack_to_df(masked_gom, var_name)
  masked_trawl <- mask_shape(masking_var, trawl_full)
  masked_trawl <- stack_to_df(masked_trawl, var_name)
  masked_dfo   <- mask_shape(masking_var, dfo_area)
  masked_dfo   <- stack_to_df(masked_dfo, var_name)
  
  # Put in list, combine
  bind_rows(list(
    "Gulf of Maine"  = masked_gom,
    "US Survey Area" = masked_trawl,
    "Canadian Survey Area"   = masked_dfo), .id = "Region")
  
}, .id = "cmip_id")





####  Check & Save  ####

# Combine and add metadata
var_combined <- bind_rows(list(masked_var, masked_percentiles)) %>% 
  mutate(
    data_source = case_when(
      str_detect(cmip_id, "percentile") ~ "Ensemble Data",
      str_detect(cmip_id, "mean") ~ "Ensemble Data",
      str_detect(cmip_id, "historic") ~ "CMIP6 Historical",
      TRUE ~ "CMIP6 SSP"),
    ensemble_statistic = case_when(
      str_detect(cmip_id, "95th") ~ "95th Percentile",
      str_detect(cmip_id, "5th") ~ "5th Percentile",
      str_detect(cmip_id, "mean") ~ "Ensemble Mean",
      TRUE ~ "Individual CMIP6 Output"),
    ensemble_statistic = factor(ensemble_statistic, 
                                levels = c("Individual CMIP6 Output",
                                           "5th Percentile",
                                           "Ensemble Mean",
                                           "95th Percentile"))
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
  facet_wrap(~Region, nrow = 1) +
  labs(x = "", y = var_label, color = "Ensemble Statistic") +
  guides(color = guide_legend(title.position = "top", title.hjust = 0.5, nrow = 2),
         alpha = "none",
         size = "none")





####  Export for later:

# Folder To Put Table
table_folder <- box_path("res", str_c("CMIP6/SSP5_85/BiasCorrected/TimeseriesData/"))

# File Name
table_name <- str_c(table_folder, "CMIP6_bias_corrected_regional_", var_select, ".csv")


# Save it
write_csv(x = var_combined, file = table_name)










