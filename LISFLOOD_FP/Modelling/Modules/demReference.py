# Prepare packages -----------------------------------------------------------------------------------------------------
from folder import *                       # For paths of sub-folders

# Packages for converting lidar point cloud las/laz file into DEM raster netcdf file
import numpy as np                         # For all calculation and data array/matrices manipulation
import os                                  # For manipulating paths (create/change/move folders) and to execute commands
import pathlib                             # For manipulating the directory path
import shapely.geometry                    # For creating a polygon object
import geopandas                           # For creating a series to store the shapely geometry object
import shutil                              # For copying file/folder
import xarray                              # For creating and manipulating an array of raster
import json                                # For creating json text data format
from geofabrics import processor           # For executing the command of converting lidar point cloud into a DEM raster
# ----------------------------------------------------------------------------------------------------------------------

def distance_calculation(point_1, point_2):
    """This function is to calculate the distance of two points using Euclidean method.

    -----------
    Reference:  https://orion.math.iastate.edu/dept/links/formulas/form2.pdf
    -----------

    -----------
    Arguments:
                point_1:
                (list)
                            First point
                point_2:
                (list)
                            Second point
    -----------

    -----------
    Returns:
                distance:
                (float)
                            Euclidean distance between two points
    -----------

    """
    # Calculating the euclidean distance between two points
    return np.sqrt(np.power((point_1[0] - point_2[0]), 2) + np.power((point_1[1] - point_2[1]), 2))


def get_divisible_number(changed_number, divisible_number=16):
    """This function is to find the padding number that is divisible by 16

    -----------
    Reference:  https://stackoverflow.com/questions/8002217/how-do-you-check-whether-a-number-is-divisible-by-another-number-python
    -----------

    -----------
    Arguments:
                changed_number:
                (float)
                                    Number that needs changing into the number that can be divisible by 16
                divisible_number:
                (float)
                                    Default is 16
    -----------

    -----------
    Return:
                changed_number:
                (float)
                                    A new number that can divisible by 16
    -----------
    """
    # Get the number that is divisible by 16
    while changed_number % divisible_number != 0:
        changed_number = changed_number + 1

    return changed_number


def padding_combination(coordinates_func, addition):
    """This function is to create the coordinates of padding/grid size

    -----------
    References:
                None.
    -----------

    -----------
    Arguments:
                coordinates_func:
                (list)
                                        A list of xmin, ymin, xmax, ymax
                addition:
                (int)
                                        A number used to extend the padding
    -----------

    -----------
    Returns:
                (list)
                                        A list of x min, x max, y min and y max
    -----------

    """
    # Get xmin, ymin, xmax, ymax
    x_min = coordinates_func[0]
    y_min = coordinates_func[1]
    x_max = coordinates_func[2]
    y_max = coordinates_func[3]

    # Calculate the center of padding
    center_x_func = ((x_min + x_max) / 2)
    center_y_func = ((y_min + y_max) / 2)
    center_padding_func = np.array([center_x_func, center_y_func])

    # Calculate the radius
    # Choose randomly a corner coordinates of rectangle point cloud to calculate the distance with center point
    radius = distance_calculation(np.array([x_min, y_min]), center_padding_func)

    # Get the coordinates of middle points of 4 sides of the rectangle point cloud
    middle_point_1 = np.array([x_min, center_y_func])
    middle_point_2 = np.array([center_x_func, y_max])
    middle_point_3 = np.array([x_max, center_y_func])
    middle_point_4 = np.array([center_x_func, y_min])

    # Calculate the distances between these points and the center point
    distance_center_1 = distance_calculation(middle_point_1, center_padding_func)
    distance_center_2 = distance_calculation(middle_point_2, center_padding_func)
    distance_center_3 = distance_calculation(middle_point_3, center_padding_func)
    distance_center_4 = distance_calculation(middle_point_4, center_padding_func)

    # Calculate the distances between the above middle points and the middle points of 4 sides of the padding
    distance_1 = radius - distance_center_1
    distance_2 = radius - distance_center_2
    distance_3 = radius - distance_center_3
    distance_4 = radius - distance_center_4

    # Convert into numbers which will be divisible by 16
    distance_padding_1 = get_divisible_number(np.ceil(distance_1))
    distance_padding_2 = get_divisible_number(np.ceil(distance_2))
    distance_padding_3 = get_divisible_number(np.ceil(distance_3))
    distance_padding_4 = get_divisible_number(np.ceil(distance_4))

    # Calculate the coordinates of the middle points of 4 sides of the padding
    padding_middle_point_1 = np.array((middle_point_1[0] - distance_padding_1, center_y_func))
    padding_middle_point_2 = np.array((center_x_func, middle_point_2[1] + distance_padding_2))
    padding_middle_point_3 = np.array((middle_point_3[0] + distance_padding_3, center_y_func))
    padding_middle_point_4 = np.array((center_x_func, middle_point_4[1] - distance_padding_4))

    # Calculate the coordinates of 4 corners of the padding
    padding_x_min = np.min(np.array([padding_middle_point_1[0],
                                     padding_middle_point_2[0],
                                     padding_middle_point_3[0],
                                     padding_middle_point_4[0]])) - addition

    padding_x_max = np.max(np.array([padding_middle_point_1[0],
                                     padding_middle_point_2[0],
                                     padding_middle_point_3[0],
                                     padding_middle_point_4[0]])) + addition

    padding_y_min = np.min(np.array([padding_middle_point_1[1],
                                     padding_middle_point_2[1],
                                     padding_middle_point_3[1],
                                     padding_middle_point_4[1]])) - addition

    padding_y_max = np.max(np.array([padding_middle_point_1[1],
                                     padding_middle_point_2[1],
                                     padding_middle_point_3[1],
                                     padding_middle_point_4[1]])) + addition

    return [padding_x_min, padding_y_min, padding_x_max, padding_y_max]


def dem_raster_reference(transformation_selection,
                         resolution_func, chunk_size_func, processor_func,
                         filename_lidar_opentopo, padding_func, lidar_dataset_name, padding=True):
    """ This function is to create a raster file from a las/laz file

    -----------
    References: https://github.com/rosepearson/GeoFabrics
    -----------

    -----------
    Arguments:
                transformation_selection:
                (string)
                                            "r" means rotation
                                            "t" means translation
                                            "c" means combination
                resolution_func:
                (int or float)
                                            resolution value in meter
                chunk_size_func:
                (int)
                                            Size value of chunk
                processor_func:
                (int)
                                            Number of processor
                filename_lidar_opentopo:
                (string)
                                            Name of the folder storing LiDAR data downloaded from Opentopography
                padding_func:
                (list)
                                            A list of x min, x max, y min and y max
                lidar_dataset_name:
                (string)
                                            LiDAR name
                padding:
                (boolean)
                                            if True, padding, else, no padding. Default is True
    -----------

    -----------
    Returns:
                None.
    -----------

    """
    # Set up the path for transformation_selection
    if transformation_selection == 'r':
        nc_raster_path_func = rotated_nc_raster_path
        elements_path_func = rotated_elements_path
        transformed = "rotated"
    elif transformation_selection == 't':
        nc_raster_path_func = translated_nc_raster_path
        elements_path_func = translated_elements_path
        transformed = "translated"
    else:
        nc_raster_path_func = combined_nc_raster_path
        elements_path_func = combined_elements_path
        transformed = "combined"

    # Name of output
    if padding:
        output_name = "generated_dem_reference.nc"
        catchment_boundary = "catchment_boundary_reference"
    else:
        output_name = "generated_dem_no_padding.nc"
        catchment_boundary = "catchment_boundary_no_padding"

    # Set crs
    h_crs = 2193
    v_crs = 7839

    # Create folder for storing data materials
    data_dir = pathlib.Path(os.getcwd()) / pathlib.Path(f"{elements_path_func}\\data_{transformed}_ref")
    if not os.path.exists(data_dir):
        os.mkdir(data_dir)

    # Bounding box - Part 1
    x0 = padding_func[0]
    y0 = padding_func[1]
    x1 = padding_func[2]
    y1 = padding_func[3]
    catchment = shapely.geometry.Polygon([(x0, y0), (x1, y0), (x1, y1), (x0, y1)])
    catchment = geopandas.GeoSeries([catchment])
    catchment = catchment.set_crs(h_crs)

    # Bounding box - Part 2
    catchment_path = pathlib.Path(f"{elements_path_func}\\data_{transformed}_ref\\{catchment_boundary}")
    catchment.to_file(catchment_path)
    shutil.make_archive(base_name=catchment_path, format='zip', root_dir=catchment_path)
    shutil.rmtree(catchment_path)

    # Reference DEM
    grid_dem_x, grid_dem_y = np.meshgrid(np.arange(np.floor(x0), np.ceil(x1), resolution_func),
                                         np.arange(np.floor(y0), np.ceil(y1), resolution_func))
    grid_dem_z = np.full(grid_dem_x.shape, fill_value=0, dtype=np.float64)
    dem = xarray.DataArray(grid_dem_z, coords={'x': grid_dem_x[0], 'y': grid_dem_y[:, 0]}, dims=['y', 'x'],
                           attrs={'scale_factor': 1.0, 'add_offset': 0.0})
    dem.rio.write_crs(h_crs)
    dem.rio.to_raster(f"{elements_path_func}\\data_{transformed}_ref\\reference_dem.nc")

    # Design JSON instructions
    instructions = {"instructions":
        {
            "output": {
                "crs": {
                    "horizontal": h_crs,
                    "vertical": v_crs
                },
                "grid_params": {
                    "resolution": resolution_func
                }
            },
            "apis": {
                "open_topography": {
                    f"{lidar_dataset_name}": {
                        "crs": {
                            "horizontal": h_crs,
                            "vertical": v_crs
                        }
                    }
                }
            },
            "processing": {
                "chunk_size": chunk_size_func,
                "number_of_cores": processor_func
            },
            "data_paths": {
                "local_cache": f"{original_lidar_path}\\{filename_lidar_opentopo}",
                "catchment_boundary": f"{elements_path_func}\\data_{transformed}_ref\\{catchment_boundary}.zip",
                "land": f"{elements_path_func}\\data_{transformed}_ref\\{catchment_boundary}.zip",
                "reference_dems": [f"{elements_path_func}\\data_{transformed}_ref\\reference_dem.nc"],
                "dense_dem_extents": f"{elements_path_func}\\data_{transformed}_ref\\extents.geojson",
                "result_dem": f"{nc_raster_path_func}\\{output_name}"
            },
            "general": {
                "set_dem_shoreline": True,
                "drop_offshore_lidar": False,
                "lidar_classifications_to_keep": [2],
                "interpolation_missing_values": True
            }
        }
    }

    # Save the instructions
    with open(f"{elements_path_func}\\data_{transformed}_ref\\instructions.json", "w") as instruction:
        json.dump(instructions, instruction)

    # Create DEM raster
    runner = processor.DemGenerator(instructions)
    runner.run()