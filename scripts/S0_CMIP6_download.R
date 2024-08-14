# Load necessary libraries
library(httr)
library(raster)
library(furrr)
library(future)
library(tidyr)  # For expand_grid

# Set up 6 parallel workers
plan(multisession, workers = 6)

# Define the base URL for WorldClim CMIP6 data
base_url <- "https://geodata.ucdavis.edu/cmip6/30s/"

# Define the list of future climate scenarios and models
scenarios <- c("ssp126", "ssp245", "ssp370", "ssp585")
models <- c("ACCESS-CM2", "CMCC-ESM2", "EC-Earth3-Veg", "FIO-ESM-2-0", "GFDL-ESM4", "GISS-E2-1-G", "HadGEM3-GC31-LL", "INM-CM5-0", "IPSL-CM6A-LR", "MIROC6", "MPI-ESM1-2-HR", "MRI-ESM2-0", "UKESM1-0-LL")

# Define the time periods
periods <- c("2021-2040", "2041-2060", "2061-2080", "2081-2100")

# Define the variable to download (e.g., temperature and precipitation)
variables <- c("bioc")
resolution <- c("30s")

# Define the base directory for organized downloads
base_dir <- "data/raw/CMIP6/"

# Function to download a single file
download_file <- function(scenario, model, period, variable, resolution) {
  # Create the directory structure if it doesn't exist
  dir_path <- paste0(base_dir, period, "/", scenario)
  if (!dir.exists(dir_path)) {
    dir.create(dir_path, recursive = TRUE)
  }
  
  # Construct the file name and URL
  file_name <- paste0("wc2.1_", resolution, "_", variable, "_", model, "_", scenario, "_", period, ".tif")
  file_url <- paste0(base_url, model, "/", scenario, "/", file_name)
  file_path <- paste0(dir_path, "/", file_name)
  
  # Check if the file already exists
  if (file.exists(file_path)) {
    message("Already exists: ", file_name, " in ", dir_path)
    return(NULL)
  }
  
  # Inform the user about the current download
  message("Downloading: ", file_name)
  
  # Download the file with error handling and progress bar
  tryCatch({
    response <- GET(file_url, write_disk(file_path, overwrite = TRUE), progress())
    if (http_status(response)$category == "Success") {
      message("Downloaded: ", file_name, " to ", dir_path)
    } else {
      message("Failed to download: ", file_name, " - HTTP status: ", http_status(response)$message)
    }
  }, error = function(e) {
    message("Failed to download: ", file_name, " - ", e$message)
  })
}

# Create a data frame with all combinations of parameters
download_tasks <- expand_grid(
  scenario = scenarios,
  model = models,
  period = periods,
  variable = variables
)

# Perform parallel downloads
future_map(seq_len(nrow(download_tasks)), function(i) {
  download_file(
    scenario = download_tasks$scenario[i],
    model = download_tasks$model[i],
    period = download_tasks$period[i],
    variable = download_tasks$variable[i],
    resolution = resolution
  )
}, .progress = TRUE)
