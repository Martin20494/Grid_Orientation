# Prepare packages -----------------------------------------------------------------------------------------------------

from folder import *                                                # For paths of sub-folders

import pandas as pd                                                 # For dataframe manipulation

import matplotlib.pyplot as plt                                     # For plotting

from randomisation import random_transformation                     # For collecting transformation values

from rasterPolygon import raster_array_to_shapefile                 # For converting points from raster to shapefile

from depthValue import run_depth_value_extraction                   # For getting depth values of all simulations

from fileWriting import csv_generation                              # For writing files into csv and raster

from checkingDifference import get_diff, plot_diff, \
                               difference_information               # For checking the differences between 0 and 90

from runStatistic import calculation_dict                           # For generating statistical dictionary

from colorSelection import *                                        # For color plotting

from statisticalPlot import plotting_map, \
                            plotting_histogram, \
                            plot_area, \
                            scatter_area_transformation             # For plotting

from savePlot import save_plot                                      # For saving plots
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

# Random transformation value
ran_trans = random_transformation(
    1,
    1,
    [0, 91, 2],
    [0, 1, 1],
    [0, 1, 1],
    'systematic',
    False
)
# ----------------------------------------------------------------------------------------------------------------------


# Clip dataset and write out point polygon -----------------------------------------------------------------------------
clipped_dataset = raster_array_to_shapefile(transform_selection, "angle_0_x_0_y_0", time_extract)
# ----------------------------------------------------------------------------------------------------------------------


# Depth values extraction and write into csv ---------------------------------------------------------------------------
# Get data dictionary with x, y coordinate values
clipped_data_coordinate = {"x_coord": clipped_dataset[:, 0],
                           "y_coord": clipped_dataset[:, 1]}

# Get depth data list of all simulations
run_depth_value_extraction(ran_trans,
                           clipped_data_coordinate,
                           transform_selection,
                           time_extract)

# Write data into pandas dataframe
full_data = pd.DataFrame(data=clipped_data_coordinate)

# Write data into csv file
csv_generation(transform_selection, full_data)
# ----------------------------------------------------------------------------------------------------------------------


# Call data from csv ---------------------------------------------------------------------------------------------------
full_data = pd.read_csv(fr"{csv_combination}\un_combined_file_1.csv")
# ----------------------------------------------------------------------------------------------------------------------


# Check differences 0 vs 90 degrees ------------------------------------------------------------------------------------
# Get a list of differences
diff_list = get_diff(full_data)

# Draw histogram
fig, ax = plt.subplots(figsize=(15, 15))
plot_diff(diff_list, ax)

# Get difference information
difference_information(
    full_data,
    diff_list
)

# ----------------------------------------------------------------------------------------------------------------------


# Statistical calculation ----------------------------------------------------------------------------------------------

# Generate calculation dictionary
statistic_dict = calculation_dict(
    transform_selection, full_data,
    resolution, filter_rate
)
# End statistical calculation ------------------------------------------------------------------------------------------



# Plot -----------------------------------------------------------------------------------------------------------------

# Set up statistical selection and color selection
statistical_selection = 'mean'
color_selection = hex_list9

# Set up plotting background -------------------------
fig, ax = plt.subplots(figsize=(20, 20))

# Plot map of statistical selection
plotting_map(
    statistic_dict[statistical_selection],
    resolution,
    statistical_selection,
    ax,
    fig,
    color_selection,
    "horizontal",
    "max"
)

# Save plot
save_plot(
    transform_selection,
    fig,
    statistical_selection,
    'png',
    50
)

# ---------------------------------------------------


# Set up plotting background ------------------------
fig, ax = plt.subplots(figsize=(20, 20))

# Plot histogram of statistical selection
plotting_histogram(
    statistic_dict[statistical_selection],
    resolution,
    statistical_selection,
    plt, ax,
    color_selection,
    0,
    11
)

# Save plot
save_plot(
    transform_selection,
    fig,
    f"{statistical_selection} histogram",
    'png',
    50
)

# ---------------------------------------------------




# Set up plotting background ------------------------
fig, ax = plt.subplots(figsize=(20, 20))

# Plot area distribution
plot_area(
    ax,
    statistic_dict['area'],
    resolution,
    'b'
)

# Save plot
save_plot(
    transform_selection,
    fig,
    f"area histogram",
    'png',
    50
)

# ---------------------------------------------------




# Set up plotting background ------------------------
fig, ax = plt.subplots(figsize=(20, 20))

# Plot area scatterplot
scatter_area_transformation(
    transform_selection,
    resolution,
    ax,
    statistic_dict['area']
)

# Save plot
save_plot(
    transform_selection,
    fig,
    f"area scatterplot",
    'png',
    50
)

# ----------------------------------------------------------------------------------------------------------------------