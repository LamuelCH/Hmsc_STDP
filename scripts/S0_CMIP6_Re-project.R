# Load necessary libraries
library(terra)
library(dplyr)

# Define directories
localDir <- "."
dataDir <- file.path(localDir, "data")
resultDir <- file.path(localDir, "results")

options(mc.cores = parallel::detectCores()-1)

# Define the list of SSPs, periods, and models
ssps <- c("ssp126", "ssp245", "ssp370", "ssp585")
periods <- c("2021-2040", "2041-2060", "2061-2080", "2081-2100")
models <- c("ACCESS-CM2", "CMCC-ESM2", "EC-Earth3-Veg", "FIO-ESM-2-0", 
            "GFDL-ESM4", "GISS-E2-1-G", "HadGEM3-GC31-LL", "INM-CM5-0", 
            "IPSL-CM6A-LR", "MIROC6", "MPI-ESM1-2-HR", "MRI-ESM2-0", "UKESM1-0-LL")

# Define function to process and save rasters
process_and_save_rasters <- function(ssp, period, models, dataDir) {
  aoi.pre <- ext(140, 150, -45, -38)
  aoi <- ext(225231.8, 629809.5, 5151752.7, 5660841.6) # Adjust as needed
  target_crs <- "EPSG:32755"  # UTM 55S
  target_resolution <- 1000  # 1km
  
  # Create an empty raster layer for projection
  raster <- rast(ext = aoi,
                 crs = target_crs,
                 resolution = target_resolution)
  
  env_tasMask <- rast("data/env_tasMask_1km.tif") %>% terra::project(., raster)
  
  # Loop through models
  for (model in models) {
    # Build file path
    file_path <- file.path(dataDir, "raw", "CMIP6", period, ssp, paste0("wc2.1_30s_bioc_", model, "_", ssp, "_", period, ".tif"))
    
    # Check if the file exists
    if (file.exists(file_path)) {
      tryCatch({
        # Load and process raster
        raster_data <- rast(file_path) %>% 
          crop(aoi.pre) %>% 
          terra::project(., env_tasMask) * env_tasMask  # Mask raster
        
        # Rename layers
        names(raster_data) <- paste0("bio", 1:19)
        
        # Define output file path
        output_file_path <- file.path("data/CMIP6/", period, ssp, paste0("wc2.1_30s_bioc_", model, "_", ssp, "_", period, "_1km.tif"))
        
        # Create directories if they don't exist
        dir.create(dirname(output_file_path), showWarnings = FALSE, recursive = TRUE)
        
        # Save the processed raster
        writeRaster(raster_data, output_file_path, overwrite = TRUE)
        
      }, error = function(e) {
        # Handle the error and continue with the next file
        warning(paste("Error processing file:", file_path, "\n", e$message))
      })
    } else {
      warning(paste("File not found:", file_path))
    }
  }
}

# Loop through SSPs and periods
for (ssp in ssps) {
  for (period in periods) {
    process_and_save_rasters(ssp, period, models, dataDir)
  }
}
