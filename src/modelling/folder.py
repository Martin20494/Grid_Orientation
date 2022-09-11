# Prepare packages --------------------------------------------------------------------------------
import pathlib                              # For creating and manipulating the directory path
import os                                   # For manipulating the directory and to execute commands
from versionModule import version          # For changing version
# -------------------------------------------------------------------------------------------


# Set up PATH -------------------------------------------------------------------------------
# Assign the path to the variable
MAIN_DIR = f"S:\\new_versions\\{version}"

# Create header path
pathlib.Path(f"{MAIN_DIR}").mkdir(parents=True, exist_ok=True)

# Change the directory
os.chdir(f"{MAIN_DIR}")
# -------------------------------------------------------------------------------------------


# Create paths #############################################################################
# ------------------------------- This is for TRANSFORMATION PART --------------------------

## General data: Folder stores all files that need to download
other_data = f"{MAIN_DIR}\\other_data"

## 0_lidar_data: Folder stores downloaded lidar point cloud data from Open Topography website
# Reference: https://portal.opentopography.org/datasets
original_lidar_path = f"{MAIN_DIR}\\0_lidar_data"

## 1_transformation: Folder stores transformed lidar point cloud data
transformed_lidar_path = f"{MAIN_DIR}\\1_transformation"

## 2_raster: Folder stores DEM rasters
# DEM -------------
# Netcdf
transformed_dem_nc_path = f"{MAIN_DIR}\\2_raster\\dem\\netcdf"
# GeoTiff
transformed_dem_tiff_path = f"{MAIN_DIR}\\2_raster\\dem\\tiff"
# Ascii
transformed_dem_asc_path = f"{MAIN_DIR}\\2_raster\\dem\\ascii"

# MANNING'S N -------
# Roughness
transformed_roughness_nc_path = f"{MAIN_DIR}\\2_raster\\n\\roughness"
# Netcdf
transformed_n_nc_path = f"{MAIN_DIR}\\2_raster\\n\\netcdf"
# GeoTiff
transformed_n_tiff_path = f"{MAIN_DIR}\\2_raster\\n\\tiff"
# Ascii
transformed_n_asc_path = f"{MAIN_DIR}\\2_raster\\n\\ascii"

# STARTDEPTH ------
# Netcdf
transformed_startdepth_nc_path = f"{MAIN_DIR}\\2_raster\\startdepth\\netcdf"
# GeoTiff
transformed_startdepth_tiff_path = f"{MAIN_DIR}\\2_raster\\startdepth\\tiff"
# Ascii
transformed_startdepth_asc_path = f"{MAIN_DIR}\\2_raster\\startdepth\\ascii"

# SHAPEFILE of PADDING
transformed_shapefile = f"{MAIN_DIR}\\2_raster\\shapefile"


## 3_LISFLOOD: Folder stores DEM rasters
# Parameters
transformed_FPpara_path = f"{MAIN_DIR}\\3_LISFLOOD_FP\\transformed_para"

# Outputs
transformed_FPoutput_path = f"{MAIN_DIR}\\3_LISFLOOD_FP\\transformed_output"


## 4_un_transformation
extracted_flowdepth = f"{MAIN_DIR}\\4_untransformation\\flowdepth"
untransformed_flowdepth = f"{MAIN_DIR}\\4_untransformation\\untransformed_flowdepth"


# ----------------------------------------- This is for ANALYSIS PART -----------------------
## 5_results
other_untransformation = f"{MAIN_DIR}\\5_analysis\\untransformed_other"

# Variation
csv_untransformation = f"{MAIN_DIR}\\5_analysis\\untransformed_csv"
raster_untransformation = f"{MAIN_DIR}\\5_analysis\\untransformed_raster"
plot_untransformation = f"{MAIN_DIR}\\5_analysis\\untransformed_plot"

# Impact
onepolygon_untransformation = f"{MAIN_DIR}\\5_analysis\\untransformed_impact\\onepolygon_nobackground"
oneraster_untransformation = f"{MAIN_DIR}\\5_analysis\\untransformed_impact\\oneraster_nobackground"
# -------------------------------------------------------------------------------------------

# Ending path creation ######################################################################


# Create a list of all sub-folders' paths --------------------------------------------------
necessary_files = [
    # Data needs downloading
    other_data,

    # 0_lidar_data
    original_lidar_path,

    # 1_transformation
    transformed_lidar_path,

    # 2_dem_raster
    # DEM
    transformed_dem_nc_path,
    transformed_dem_tiff_path,
    transformed_dem_asc_path,
    # MANNING'S N
    transformed_roughness_nc_path,
    transformed_n_nc_path,
    transformed_n_tiff_path,
    transformed_n_asc_path,
    # STARTDEPTH
    transformed_startdepth_nc_path,
    transformed_startdepth_tiff_path,
    transformed_startdepth_asc_path,
    # SHAPEFILE
    transformed_shapefile,

    # 3_LISFLOOD
    transformed_FPpara_path,
    transformed_FPoutput_path,

    # 4_untransformation
    extracted_flowdepth,
    untransformed_flowdepth,

    # 5_results
    other_untransformation,
    # Variation
    csv_untransformation,
    raster_untransformation,
    plot_untransformation,
    # Impact
    onepolygon_untransformation,
    oneraster_untransformation
]
# -------------------------------------------------------------------------------------------


# Execute all sub-folders -------------------------------------------------------------------
for each_folder in necessary_files:
    pathlib.Path(f"{each_folder}").mkdir(parents=True, exist_ok=True)
# -------------------------------------------------------------------------------------------