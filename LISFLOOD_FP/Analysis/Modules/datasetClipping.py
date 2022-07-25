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

# For polygon clipping
import geopandas as gpd                            # For creating polygon
# ----------------------------------------------------------------------------------------------------------------------


# Remove all warnings  -------------------------------------------------------------------------------------------------
# Reference: https://numba.pydata.org/numba-doc/dev/reference/deprecation.html#suppressing-deprecation-warnings
from numba.core.errors import NumbaDeprecationWarning, NumbaPendingDeprecationWarning, NumbaWarning
import warnings

warnings.simplefilter('ignore', category=NumbaDeprecationWarning)
warnings.simplefilter('ignore', category=NumbaPendingDeprecationWarning)
warnings.simplefilter('ignore', category=NumbaWarning)
# ----------------------------------------------------------------------------------------------------------------------


def array_creation(data_array, value, switch=False):
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
                switch:
                (boolean)
                                        To switch the row and column of array
                                        True is switching
                                        False is not switching (default)
    -----------

    -----------
    Returns:
                new_array:
                (array)
                                        New clipped array
    -----------

    """
    # Read x, y values into arrays
    arr_x = data_array.x.values
    arr_y = data_array.y.values

    if switch:
        # Create zero array
        new_array = np.zeros((arr_y.shape[0], arr_x.shape[0]))
    else:
        # Create zero array
        new_array = np.zeros((arr_x.shape[0], arr_y.shape[0]))

    # Get full number of x or y values
    if value == "x":
        for i in range(arr_x.shape[0]):
            for j in range(arr_y.shape[0]):
                new_array[j, i] = arr_x[i]
        return new_array

    else:

        for i in range(arr_x.shape[0]):
            for j in range(arr_y.shape[0]):
                new_array[j, i] = arr_y[j]
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
                                            Amount of time that flood model predicted
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
    array_x = array_creation(dataset_rioxarray, 'x', False)
    array_y = array_creation(dataset_rioxarray, 'y', False)
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


def clip(transformation_selection, dataset, adjusted_value_list=[0, 0, 0, 0], clip_using_shape_file=False):
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
                adjusted_value_list:
                (list)
                                            List contain values to change boundaries
                clip_using_shape_file:
                                            Selection for clipping by using shapfile or manually
                                            True means using shapefile and clipping out all -999/nodata values
                                            False means clipping manually (following rectangle shape)
                                            Default is False
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

    # Get new xmin, xmax, ymin, ymax
    new_xmin = x_min + adjusted_value_list[0]
    new_ymin = y_min + adjusted_value_list[1]
    new_xmax = x_max + adjusted_value_list[2]
    new_ymax = y_max + adjusted_value_list[3]

    # Calculate coordinates of center point
    center_point = center_calculation('c', True)
    center_x = center_point[0]  # Extract x coordinate of center point
    center_y = center_point[1]  # Extract y coordinate of center point

    # Remove padding by filtering
    if clip_using_shape_file:
        # if clipping using in-rectangle shapefile
        dataset_clip = dataset[(dataset[:, 2] > -999)]
    else:
        # if clipping using rectangle shape
        dataset_clip = dataset[(new_xmin <= dataset[:, 0])
                               & (new_xmax >= dataset[:, 0])
                               & (new_ymin <= dataset[:, 1])
                               & (new_ymax >= dataset[:, 1])]

    # Adjust x and y to get flowdepth values later (to avoid the case one point with two flowdepth values)
    # Rotate -0.000001 degree
    adjusted_dataset_clip1 = wrapping_point_rotation(dataset_clip, -0.000001, center_x, center_y, 0)
    # Translate -0.000001 meter
    adjusted_dataset_clip2 = adjusted_dataset_clip1.copy()
    adjusted_dataset_clip2[:, 0] = adjusted_dataset_clip1[:, 0] - 0.000001
    adjusted_dataset_clip2[:, 1] = adjusted_dataset_clip1[:, 1] - 0.000001
    adjusted_dataset_clip2[:, 2] = adjusted_dataset_clip1[:, 2]

    return adjusted_dataset_clip2


def clipping_using_shapefile(transformation_selection, polygon_coord,
                            number_simulation, time_extract_func,
                            value_clipping):
    """This function is to clip the raster using created shapefile

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
                polygon_coord:
                (shapely polygon format)
                                            None means polygon shapfile was created manually (e.g. using QGIS)
                                            or shapely polygon format (e.g. Polygon([(0, 1), (1, 1), (1, 0), (0, 0)]))
                number_simulation:
                (int)
                                            Ordinal number of simulation
                time_extract_func:
                (int)
                                            Amount of time that flood model predicte
                value_clipping:
                (float or int)
                                            This value is used for clipping, in which all values within polygon will
                                            be changed into this value (e.g. -999)
    -----------

    -----------
    Returns:
                None.
    -----------

    """
    # Set up the path for transformation_selection
    if transformation_selection == 'r':
        clipping_shapefile = polygon_rotation
        transformed = "rotated"
        flowdepth_path = rotated_flowdepth
    elif transformation_selection == 't':
        clipping_shapefile = polygon_translation
        transformed = "translated"
        flowdepth_path = translated_flowdepth
    else:
        clipping_shapefile = polygon_combination
        transformed = "combined"
        flowdepth_path = combined_flowdepth

    # Create polygon based on the given list of polygon coordinates
    if polygon_coord is not None:
        polygon_data = gpd.GeoDataFrame(geometry=[polygon_coord], crs=2193)
        polygon_data.to_file(fr"{clipping_shapefile}\\polygon_clipping.shp")
    else:
        pass

    # Get the path of polygon_clipping.shp
    polygon_clipping_path = fr"{clipping_shapefile}\\polygon_clipping.shp"

    # Convert raster from netcdf to geotiff
    # Get raster from 5_flowdepth file
    path = fr"{flowdepth_path}\\flowdepth_{transformed}_{number_simulation}_at_{time_extract_func}.nc"
    raster_nc = rioxarray.open_rasterio(path)
    # Convert to geotiff file and save to polygon file in 7_results
    raster_nc.rio.to_raster(
        fr"{clipping_shapefile}\\converted_flowdepth_{transformed}_{number_simulation}_at_{time_extract_func}.tiff"
    )
    # Get raster geotiff path
    raster_tiff_path = fr"{clipping_shapefile}\\converted_flowdepth_{transformed}_{number_simulation}_at_{time_extract_func}.tiff"

    # Clip the raster by changing all values into a specific value (e.g. -999)
    # Build up command
    command = fr"gdal_rasterize -i -burn -999 {polygon_clipping_path} {raster_tiff_path}"
    # Execute the command
    os.system(command)

    # Convert value_clipping to nodata value
    # Read the geotiff file
    raster_tiff = rioxarray.open_rasterio(
        raster_tiff_path
    )
    # Convert
    nodata_raster_tiff = raster_tiff.rio.write_nodata(value_clipping, inplace=True)
    # Add crs
    crs_nodata_raster_tiff = nodata_raster_tiff.rio.write_crs(2193)
    # Write to file
    crs_nodata_raster_tiff.rio.to_raster(fr"{flowdepth_path}\\flowdepth_{transformed}_{number_simulation}_clipped_at_{time_extract_func}.nc")

