# Prepare packages -----------------------------------------------------------------------------------------------------
import pathlib  # For creating and manipulating the directory path
import os  # For manipulating the directory and to execute commands
# ----------------------------------------------------------------------------------------------------------------------


# Set up 'header' path -------------------------------------------------------------------------------------------------
drive = "S"
main_folder = "LISFLOOD"
version = "version_25"
header = fr"{drive}:\\{main_folder}\\{version}"

# Create header path
pathlib.Path(f"{header}").mkdir(parents=True, exist_ok=True)

# Change the directory
os.chdir(f"{header}")
# ----------------------------------------------------------------------------------------------------------------------


# Create paths  ########################################################################################################
# -------------------------------------- This is for TRANSFORMATION PART -------------------------------------------

# 1_open_topography: Folder stores ownloaded lidar point cloud data from open topography website
# Reference: https://portal.opentopography.org/datasets
original_lidar_path = f"{header}\\0_open_topography"

# 2_transformation: Folder stores transformed lidar point cloud data
rotated_lidar_path = f"{header}\\2_transformation\\rotated_lidar"
translated_lidar_path = f"{header}\\2_transformation\\translated_lidar"
combined_lidar_path = f"{header}\\2_transformation\\combined_lidar"

# 3_dem_raster: Folder stores DEM rasters
# Rotation
# Netcdf
rotated_nc_raster_path = f"{header}\\3_dem_raster\\rotated_dem_raster\\rot_dem_raster\\rot_netcdf"
# GeoTiff
rotated_tiff_raster_path = f"{header}\\3_dem_raster\\rotated_dem_raster\\rot_dem_raster\\rot_tiff"
# Ascii
rotated_asc_raster_path = f"{header}\\3_dem_raster\\rotated_dem_raster\\rot_dem_raster\\rot_ascii"
rotated_elements_path = f"{header}\\3_dem_raster\\rotated_dem_raster\\rot_elements"

# Translation
# Netcdf
translated_nc_raster_path = f"{header}\\3_dem_raster\\translated_dem_raster\\tra_dem_raster\\tra_netcdf"
# GeoTiff
translated_tiff_raster_path = f"{header}\\3_dem_raster\\translated_dem_raster\\tra_dem_raster\\tra_tiff"
# Ascii
translated_asc_raster_path = f"{header}\\3_dem_raster\\translated_dem_raster\\tra_dem_raster\\tra_ascii"
translated_elements_path = f"{header}\\3_dem_raster\\translated_dem_raster\\tra_elements"

# Combination
# Netcdf
combined_nc_raster_path = f"{header}\\3_dem_raster\\combined_dem_raster\\com_dem_raster\\com_netcdf"
# GeoTiff
combined_tiff_raster_path = f"{header}\\3_dem_raster\\combined_dem_raster\\com_dem_raster\\com_tiff"
# Ascii
combined_asc_raster_path = f"{header}\\3_dem_raster\\combined_dem_raster\\com_dem_raster\\com_ascii"
combined_elements_path = f"{header}\\3_dem_raster\\combined_dem_raster\\com_elements"

# 4_LISFLOOD: Folder stores LISFLOOD-FP output
# LISFLOOD_FP output
rotated_FPoutput_path = f"{header}\\4_LISFLOOD_FP\\rotated_output"
translated_FPoutput_path = f"{header}\\4_LISFLOOD_FP\\translated_output"
combined_FPoutput_path = f"{header}\\4_LISFLOOD_FP\\combined_output"

# 5_flowdepth: Folder stores flowdepth extracted from LISFLOOD-FP output
rotated_flowdepth = f"{header}\\5_flowdepth\\rotated_flowdepth"
translated_flowdepth = f"{header}\\5_flowdepth\\translated_flowdepth"
combined_flowdepth = f"{header}\\5_flowdepth\\combined_flowdepth"

# 6_un_transformation: Folder stores un-transformed raster
unrotated_path = f"{header}\\6_un_transformation\\un_rotation"
untranslated_path = f"{header}\\6_un_transformation\\un_translation"
uncombined_path = f"{header}\\6_un_transformation\\un_combination"

# ----------------------------------------- This is for ANALYSIS PART ----------------------------------------------

# 7_results: Folder stores mean and standard deviation
# Un_rotation
csv_rotation = f"{header}\\7_results\\un_rotation\\unrot_csv"
raster_rotation = f"{header}\\7_results\\un_rotation\\unrot_raster"
polygon_rotation = f"{header}\\7_results\\un_rotation\\unrot_polygon"
plot_rotation = f"{header}\\7_results\\un_rotation\\unrot_plot"

# Un_translation
csv_translation = f"{header}\\7_results\\un_translation\\untra_csv"
raster_translation = f"{header}\\7_results\\un_translation\\untra_raster"
polygon_translation = f"{header}\\7_results\\un_translation\\untra_polygon"
plot_translation = f"{header}\\7_results\\un_translation\\untra_plot"

# Un_combination
csv_combination = f"{header}\\7_results\\un_combination\\uncom_csv"
raster_combination = f"{header}\\7_results\\un_combination\\uncom_raster"
polygon_combination = f"{header}\\7_results\\un_combination\\uncom_polygon"
plot_combination = f"{header}\\7_results\\un_combination\\uncom_plot"

# Ending path creation #################################################################################################


# Create a list of all sub-folders' paths ------------------------------------------------------------------------------
necessary_files = [
    # 1_open_topography
    original_lidar_path,

    # 2_transformation
    rotated_lidar_path,
    translated_lidar_path,
    combined_lidar_path,

    # 3_dem_raster
    # Rotation
    rotated_nc_raster_path,
    rotated_tiff_raster_path,
    rotated_asc_raster_path,
    rotated_elements_path,

    # Translation
    translated_nc_raster_path,
    translated_tiff_raster_path,
    translated_asc_raster_path,
    translated_elements_path,

    # Combination
    combined_nc_raster_path,
    combined_tiff_raster_path,
    combined_asc_raster_path,
    combined_elements_path,

    # 4_LISFLOOD-FP
    rotated_FPoutput_path,
    translated_FPoutput_path,
    combined_FPoutput_path,

    # 5_flowdepth
    rotated_flowdepth,
    translated_flowdepth,
    combined_flowdepth,

    # 6_un_transformation
    unrotated_path,
    untranslated_path,
    uncombined_path,

    # 7_results
    # Unrotation
    csv_rotation,
    raster_rotation,
    polygon_rotation,
    plot_rotation,

    # Untranslation
    csv_translation,
    raster_translation,
    polygon_translation,
    plot_translation,

    # Uncombination
    csv_combination,
    raster_combination,
    polygon_combination,
    plot_combination
]
# ----------------------------------------------------------------------------------------------------------------------


# Create all sub-folders -----------------------------------------------------------------------------------------------
for each_folder in necessary_files:
    pathlib.Path(f"{each_folder}").mkdir(parents=True, exist_ok=True)
# ----------------------------------------------------------------------------------------------------------------------
