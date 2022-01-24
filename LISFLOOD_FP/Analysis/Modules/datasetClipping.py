# Prepare packages -----------------------------------------------------------------------------------------------------
from folder import *                               # For paths of sub-folders

# Packages for transformation
from numba import guvectorize, float64             # For speeding up the time of the code running

# Packages for unrotating and untranslating
from osgeo import gdal                             # For manipulating rasters (calculating centers)

# For handling raster
import numpy as np                                 # For all calculation and data array/matrices manipulation
import rioxarray                                   # For manipulating pixel values, spatial attributes, and raster files
import rasterio                                    # For manipulating raster files (mainly crs and transformation)
# ----------------------------------------------------------------------------------------------------------------------


# Remove all warnings  -------------------------------------------------------------------------------------------------
# Reference: https://numba.pydata.org/numba-doc/dev/reference/deprecation.html#suppressing-deprecation-warnings
from numba.core.errors import NumbaDeprecationWarning, NumbaPendingDeprecationWarning, NumbaWarning
import warnings

warnings.simplefilter('ignore', category=NumbaDeprecationWarning)
warnings.simplefilter('ignore', category=NumbaPendingDeprecationWarning)
warnings.simplefilter('ignore', category=NumbaWarning)
# ----------------------------------------------------------------------------------------------------------------------


def array_creation(data_array, value, shape_x, shape_y):
    """This function is to create an array of a coordinate

    -----------
    Reference:
                None.
    -----------

    -----------
    Arguments:
                data_array:
                (array)
                                        Original array with full data
                value:
                (string)
                                        The value of original array need extracting.
                                        Two main choices - 'x' or 'y'
                shape_x:
                (int)
                                        Shape of x array
                shape_y:
                (int)
                                        Shape of y array
    -----------

    -----------
    Returns:
                new_array:
                (array)
                                        New clipped array
    -----------

    """
    # Create zero array
    new_array = np.zeros((shape_x, shape_y))

    # Read x, y values into arrays
    arr_x = data_array.x.values
    arr_y = data_array.y.values

    # Get full number of x or y values
    for i in range(new_array.shape[0]):
        for j in range(new_array.shape[0]):
            if value == 'x':
                new_array[j, i] = arr_x[i]
            else:
                new_array[i, j] = arr_y[i]

    return new_array


def xyz_array(transformation_selection, number_simulation, time_extract_func):
    """This function is to create an array from flowdepth raster including x, y, z values

    -----------
    Reference:
                None.
    -----------

    -----------
    Arguments:
                transformation_selection:
                (string)
                                            "r" means rotation
                                            "t" means translation
                                            "c" means combination
                number_simulation:
                (int)
                                            Ordinal number of simulation
                time_extract_func:
                (int)
                                            Amount of time that BG_flood model predicted
    -----------

    -----------
    Returns:
                full_dataset:
                (array)
                                            A dataset contains 1D x, y, z array including padding
    -----------

    """
    # Set up the path for transformation_selection
    if transformation_selection == 'r':
        transformed = "rotated"
        flowdepth_path = rotated_flowdepth
    elif transformation_selection == 't':
        transformed = "translated"
        flowdepth_path = translated_flowdepth
    else:
        transformed = "combined"
        flowdepth_path = combined_flowdepth

    # Read the raster
    path = fr"{flowdepth_path}\\flowdepth_{transformed}_{number_simulation}_at_{time_extract_func}.nc"
    dataset_rioxarray = rioxarray.open_rasterio(path)

    # Create full number of values of x, y, z coordinates
    array_x = array_creation(dataset_rioxarray, 'x', dataset_rioxarray.shape[1], dataset_rioxarray.shape[2])
    array_y = array_creation(dataset_rioxarray, 'y', dataset_rioxarray.shape[1], dataset_rioxarray.shape[2])
    array_z = dataset_rioxarray.isel(band=0).values

    # Flatten x, y, z arrays
    flatten_x = array_x.flatten()
    print(flatten_x)
    flatten_y = array_y.flatten()
    print(flatten_y)
    flatten_z = array_z.flatten()
    print(flatten_z)

    # Put all x, y, z into one array
    full_dataset = np.vstack((flatten_x, flatten_y, flatten_z)).transpose()

    return full_dataset


# ROTATION #############################################################################################################
def center_calculation(transformation_selection, lidar=True):
    """This function is to calculate the center of reference DEM raster which will be used in the whole process

    -----------
    References: https://gis.stackexchange.com/questions/104362/how-to-get-extent-out-of-geotiff
                https://rasterio.readthedocs.io/en/latest/quickstart.html
    -----------

    -----------
    Arguments:
                transformation_selection:
                (string)
                                            "r" means rotation
                                            "t" means translation
                                            "c" means combination
                lidar:
                (boolean)
                                            True means doing rotation on LiDAR data
                                            False means doing unrotation on model outputs
    -----------

    -----------
    Return:
                (tuple)
                                            A tuple of x, y coordinates of center raster
    -----------

    """
    # Set up the path for transformation_selection
    if transformation_selection == 'r':
        transformed_nc_raster_path = rotated_nc_raster_path
    elif transformation_selection == 't':
        transformed_nc_raster_path = translated_nc_raster_path
    else:
        transformed_nc_raster_path = combined_nc_raster_path

    # Read the file using gdal
    raster_reference_func = gdal.Open(fr"{transformed_nc_raster_path}\\generated_dem_reference.nc")

    if lidar:
        # Get the center coordinates for LiDAR
        xmin, xpixel, _, ymax, _, ypixel = raster_reference_func.GetGeoTransform()
        width, height = raster_reference_func.RasterXSize, raster_reference_func.RasterYSize
        xmax = xmin + width * xpixel
        ymin = ymax + height * ypixel
        center_x_func = ((xmin + xmax) / 2)
        center_y_func = ((ymin + ymax) / 2)
        center_func = np.array([center_x_func, center_y_func])

    else:
        # Get the center coordinates for DEMs
        center_func = (raster_reference_func.RasterXSize / 2, raster_reference_func.RasterYSize / 2)
    return center_func



# Create a command for guvectorize() decorator
# In that, float64[:,:] stands for the format of the array and float64 stands for the format of other parameters
# (m,n) signature stands for the matrix of the array and () signature stands for other parameters
gu_rotation = guvectorize([(float64[:, :], float64, float64, float64, float64, float64[:, :])],
                          '(m,n),(),(),(),()->(m,n)')

def point_rotation(coordinates_func, angle, center_x_func, center_y_func, clockwise, new_coordinates_func):
    """This function is to calculate the rotated coordinates of LiDAR data

    -----------
    References: https://www.youtube.com/watch?v=RqZH-7hlI48
                https://stackoverflow.com/questions/14607640/rotating-a-vector-in-3d-space
                https://stackoverflow.com/questions/5954603/transposing-a-1d-numpy-array
                https://en.wikipedia.org/wiki/Rotation_matrix
                https://math.stackexchange.com/questions/270194/how-to-find-the-vertices-angle-after-rotation
                https://github.com/numba/numba/issues/3312
                https://numba.pydata.org/numba-doc/latest/user/vectorize.html
                https://numba.pydata.org/numba-doc/latest/cuda/ufunc.html
                http://numba.pydata.org/numba-doc/0.20.0/reference/compilation.html
                http://numba.pydata.org/numba-doc/0.12/tutorial_numpy_and_numba.html
                https://numba.pydata.org/numba-doc/dev/reference/types.html
    -----------

    -----------
    Arguments:
                coordinates_func:
                (array)
                                        An array of the coordinates of the point will be rotated
                angle:
                (int or float)
                                        The value of angle to rotate
                center_x_func:
                (float)
                                        Coordinate value of x center.
                                        The x center here was used from the center of reference DEM without padding
                center_y_func:
                (float)
                                        Coordinate value of x center.
                                        The x center here was used from the center of reference DEM without padding
                clockwise:
                (string)
                                        Rotating the points in clockwise (1) or anti_clockwise (0) directions
                new_coordinates_func:
                (array)
                                        A new array of rotated x, y, z coordinates values of a point
    -----------

    -----------
    Returns:
                (array)
                                        The function will return itself without the need of 'return' command
    -----------

    """
    # Convert degree to radian and calculate cosine and sine
    radian = np.deg2rad(angle)
    cosine = np.cos(radian)
    sine = np.sin(radian)

    # Create a for loop to manipulate each row of the array
    for i in range(coordinates_func.shape[0]):
        # Do substraction with center coordinates
        diff_x = coordinates_func[i, 0] - center_x_func
        diff_y = coordinates_func[i, 1] - center_y_func

        # Calculate the rotated point coordinates
        if clockwise == 0:  # Rotating in anti-clockwise direction
            new_coordinates_func[i, 0] = diff_x * cosine - diff_y * sine + center_x_func
            new_coordinates_func[i, 1] = diff_x * sine + diff_y * cosine + center_y_func
            new_coordinates_func[i, 2] = coordinates_func[i, 2]
        else:  # Rotating in clockwise direction
            new_coordinates_func[i, 0] = diff_x * cosine + diff_y * sine + center_x_func
            new_coordinates_func[i, 1] = diff_x * (-sine) + diff_y * cosine + center_y_func
            new_coordinates_func[i, 2] = coordinates_func[i, 2]

# Wrapping function to map later
wrapping_point_rotation = gu_rotation(point_rotation)

# END ROTATION #########################################################################################################


def clip(transformation_selection, dataset):
    """This function is to clip rasters (remove padding box)

    -----------
    References:
                None.
    -----------

    -----------
    Arguments:
                transformation_selection:
                (string)
                                            "r" means rotation
                                            "t" means translation
                                            "c" means combination
                dataset:
                (array)
                                            An array with x, y, z values
    -----------

    -----------
    Returns:
                adjusted_dataset_clip:
                                            A dataset without padding.
                                            The coordinates are adjusted 0.000001 degree and meter for value extraction
    -----------

    """
    # Set up the path for transformation_selection
    if transformation_selection == 'r':
        nc_raster_path_func = rotated_nc_raster_path
    elif transformation_selection == 't':
        nc_raster_path_func = translated_nc_raster_path
    else:
        nc_raster_path_func = combined_nc_raster_path

    # Get boundary coordinates
    raster_origin = rasterio.open(fr"{nc_raster_path_func}\\generated_dem_no_padding.nc")
    x_min, y_min, x_max, y_max = raster_origin.bounds

    # Calculate coordinates of center point
    center_point = center_calculation('c', True)
    center_x = center_point[0]  # Extract x coordinate of center point
    center_y = center_point[1]  # Extract y coordinate of center point

    # Remove padding by filtering
    dataset_clip = dataset[(x_min <= dataset[:, 0])
                           & (x_max >= dataset[:, 0])
                           & (y_min <= dataset[:, 1])
                           & (y_max >= dataset[:, 1])]

    # Adjust x and y to get flowdepth values later (to avoid the case one point with two flowdepth values)
    # Rotate -0.000001 degree
    adjusted_dataset_clip1 = wrapping_point_rotation(dataset_clip, -0.000001, center_x, center_y, 0)
    # Translate -0.000001 meter
    adjusted_dataset_clip2 = adjusted_dataset_clip1.copy()
    adjusted_dataset_clip2[:, 0] = adjusted_dataset_clip1[:, 0] - 0.000001
    adjusted_dataset_clip2[:, 1] = adjusted_dataset_clip1[:, 1] - 0.000001
    adjusted_dataset_clip2[:, 2] = adjusted_dataset_clip1[:, 2]

    return adjusted_dataset_clip2