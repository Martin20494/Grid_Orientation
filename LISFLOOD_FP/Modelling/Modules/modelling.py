# Prepare packages -----------------------------------------------------------------------------------------------------
from randomisation import random_transformation     # For collecting transformation values
from lidar import download_lidar                    # For downloading LiDAR data
from transformation import get_url                  # For getting url file for transformation command
from demReference import padding_combination, \
                         dem_raster_reference       # For padding calculation and dem reference generation
from transformation import center_calculation       # For center calculation
from execution import *                             # For executing transformed and un-transformed functions
# ----------------------------------------------------------------------------------------------------------------------


# Main parameters - Part 1----------------------------------------------------------------------------------------------
# LiDAR number file
num_lidar_file = 0

# Transformation selection
transform_selection = "c"

# Construct a list of boundary
resolution = 10
number_pixel_x = 16 * 20  # multiple of 16
number_pixel_y = 16 * 10  # multiple of 16

xmin = 1769850
ymin = 5471830
xmax = xmin + resolution * number_pixel_x
ymax = ymin + resolution * number_pixel_y

# Boundary for the area of interest
boundary_1 = [xmin, ymin, xmax, ymax]

# Padding boundary for DEMs (use this boundary to create padding)
addition_2 = 16 * 3
boundary_2 = padding_combination(boundary_1, addition_2)

# Boundary for tiles (use this boundary to download LiDAR)
addition_3 = 16 * 3
boundary_3 = padding_combination(boundary_2, addition_3)

# Chunks and processors information
size_of_chunk = 100
size_of_processor = 4

# Extracting flowdepth rate (for flowdepth_extraction() function)
flowdepth_rate = 0

# Time to extract flowdepth (for flowdepth_extraction() function)
time_extract = 12
# ----------------------------------------------------------------------------------------------------------------------


# Main parameters - Part 2----------------------------------------------------------------------------------------------
# Download LiDAR
download_lidar(boundary_3, num_lidar_file)

# Reference DEM without padding
dem_raster_reference(transform_selection,
                     resolution, size_of_chunk, size_of_processor,
                     f"lidar_{num_lidar_file}", boundary_1, False)

# Reference DEM with padding
dem_raster_reference(transform_selection,
                     resolution, size_of_chunk, size_of_processor,
                     f"lidar_{num_lidar_file}", boundary_2)

# Calculate coordinates of center point
center_point = center_calculation(transform_selection, True)
center_x = center_point[0]                     # Extract x coordinate of center point
center_y = center_point[1]                     # Extract y coordinate of center point

# Transformation values
ran_trans = random_transformation(
    1,
    1,
    [0, 1, 1],
    [0, 1, 1],
    [0, 1, 1],
    'systematic',
    False
)

# Get url list file
url_list_file = get_url(num_lidar_file)
# ----------------------------------------------------------------------------------------------------------------------

# Transform tiles
run_simulation(ran_trans,
               transform_selection, num_lidar_file,
               center_x, center_y,
               url_list_file)

# DEM generation
run_dem(ran_trans, transform_selection,
        resolution, size_of_chunk, size_of_processor,
        boundary_2)

# Values change
run_changing_value(ran_trans, transform_selection, center_x, center_y)

# Flood model
run_flood_model(ran_trans, transform_selection, resolution, center_x, center_y)

# Flowdepth extraction
run_flowdepth_extraction(ran_trans, transform_selection, time_extract)

# Untransformation
run_untransformation(ran_trans, transform_selection, time_extract)
