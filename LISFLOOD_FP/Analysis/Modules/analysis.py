# Prepare packages -----------------------------------------------------------------------------------------------------
import pandas as pd                                                 # For dataframe manipulation

import matplotlib.pyplot as plt                                     # For plotting

from folder import *                                                # For paths of sub-folders

from randomisation import random_transformation                     # For collecting transformation values

from rasterPolygon import raster_array_to_shapefile                 # For converting points from raster to shapefile

from depthValue import run_depth_value_extraction                   # For getting depth values of all simulations

from fileWriting import csv_generation, raster_generation           # For writing files into csv and raster

from checkingDiference import get_diff, plot_diff                   # For checking the differences between 0 and 90

from runStatistic import calculation_dict                           # For generating statistical dictionary

# ----------------------------------------------------------------------------------------------------------------------


# Main parameters ------------------------------------------------------------------------------------------------------
# Transformation selection
transform_selection = "c"

# Time extract
time_extract = 12

# Resolution
resolution = 10

# Filter rate
filter_rate = 0.1

# Clip dataset and write out point polygon
clipped_dataset = raster_array_to_shapefile(transform_selection, "angle_0_x_0_y_0", time_extract)

# Random transformation value
ran_trans = random_transformation(
    1,
    1,
    [0, 91, 2],
    [0, 2, 1],
    [0, 2, 1],
    'systematic',
    False
)
# ----------------------------------------------------------------------------------------------------------------------



# Depth values extraction and write into csv ---------------------------------------------------------------------------
# Get data dictionary with x, y coordinate values
clipped_data_coordinate = {"x_coord": clipped_dataset[:, 0],
                           "y_coord": clipped_dataset[:, 1]}

# Get depth data list of all simulations
data = run_depth_value_extraction(ran_trans,
                                  clipped_data_coordinate,
                                  transform_selection,
                                  time_extract)

# Write data into pandas dataframe
full_data = pd.DataFrame(data=data)

# Write data into csv file
csv_generation(transform_selection, full_data)
# ----------------------------------------------------------------------------------------------------------------------



# Check differences 0 vs 90 degrees ------------------------------------------------------------------------------------
# Get a list of differences
diff_list = get_diff(full_data)

# Draw histogram
fig, ax = plt.subplots(figsize=(15, 15))
plot_diff(diff_list, ax)
# ----------------------------------------------------------------------------------------------------------------------


# Statistical calculation ----------------------------------------------------------------------------------------------
# In case the data need loading again
# full_data = pd.read_csv(fr"{csv_combination}\un_combined_file_1.csv")

# Generate calculation dictionary
statistic_dict = calculation_dict(
    transform_selection, full_data,
    resolution, filter_rate
)
# End statistical calculation ------------------------------------------------------------------------------------------



