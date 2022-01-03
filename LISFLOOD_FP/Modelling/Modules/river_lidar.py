# Prepare packages -----------------------------------------------------------------------------------------------------

import numpy as np

import geopandas as gpd                             # For polygon manipulation
import pdal                                         # For LiDAR manipulation
import json                                         # For instruction manipulation
from scipy import interpolate                       # For interpolation

import rioxarray                                    # For raster manipulation
import xarray                                       # For raster manipulation

from transformation import read_las_file
# ----------------------------------------------------------------------------------------------------------------------


# FUNCTION TO CROP -----------------------------------------------------------------------------------------------------
def crop_lidar(shape_file, file_to_crop, file_being_cropped, inside=True):
    """This function is to crop LiDAR based on polygon. The results could be LiDAR outside or inside polygon

    -----------
    References:
                None.
    -----------

    -----------
    Arguments:
                shape_file:
                (string)
                                        The directory to polygon file
                file_to_crop:
                (string)
                                        The directory to the file that will be cropped
                file_being_cropped:
                (string)
                                        The directory to the cropped file
                inside:
                (boolean)
                                        True - remove all LiDAR outside
                                        False - remove all LiDAR inside
    -----------

    -----------
    Returns:
                None.
    -----------

    """
    # Get polygon figures
    polygon = gpd.read_file(shape_file)
    polygon_wkt = polygon['geometry'].to_wkt()

    # Develop cropping instruction
    if inside:
        instruction_crop = [
            file_to_crop,
            {
                "type": "filters.crop",
                "polygon": polygon_wkt[0],
                "distance": 0,
                "outside": False
            },
            {
                "type": "writers.las",
                "filename": file_being_cropped,
                "compression": "laszip"
            }
        ]
    else:
        instruction_crop = [
            file_to_crop,
            {
                "type": "filters.crop",
                "polygon": polygon_wkt[0],
                "distance": 0,
                "outside": True
            },
            {
                "type": "writers.las",
                "filename": file_being_cropped,
                "compression": "laszip"
            }
        ]

    pipeline_crop = pdal.Pipeline(json.dumps(instruction_crop))
    pipeline_crop.execute()

# END FUNCTION TO CROP -------------------------------------------------------------------------------------------------


# WRITE INTO LAZ/LAS FILE ----------------------------------------------------------------------------------------------

def writing_las_file(data_array_func):
    """This function is to create a new lidar las/laz file
    the same as the old lidar las/laz file with different coordinates

    -----------
    References: https://pdal.io/tutorial/las.html
    -----------

    -----------
    Arguments:
                data_array_func:
                (array)
                                            Array of raster

    -----------

    -----------
    Returns:
                None.
    -----------
    """
    # Set up horizontal and vertical crs
    h_crs = 2193
    v_crs = 7839

    # Get the shape of transformed LiDAR array
    lidar_arr_shape = data_array_func.shape[0]

    # Create empty transformed LiDAR array
    lidar_arr = np.empty([lidar_arr_shape], dtype=[('X', '<f8'), ('Y', '<f8'), ('Z', '<f8'), ('Classification', 'u1')])

    # Overwrite the X, Y, and Z values in the points array with the transformed values
    lidar_arr['X'] = (data_array_func[:, 0]).flatten()
    lidar_arr['Y'] = (data_array_func[:, 1]).flatten()
    lidar_arr['Z'] = (data_array_func[:, 2]).flatten()
    lidar_arr['Classification'] = np.array([2] * lidar_arr_shape)

    # Save the file
    pdal_pipeline_instructions = [
        {
            "type": "writers.las",
            "scale_x": 0.0001,
            "scale_y": 0.0001,
            "offset_x": 1700000,
            "offset_y": 5400000,
            "a_srs": f"EPSG:{h_crs}+{v_crs}",
            "filename": r"S:\Lidar_example\lidar_1m.laz",
            "compression": "laszip"
        }
    ]

    # Run
    pdal_pipeline = pdal.Pipeline(json.dumps(pdal_pipeline_instructions), [lidar_arr])
    pdal_pipeline.execute()

# END WRITING INTO LAZ/LAS FILE ----------------------------------------------------------------------------------------


# FUNCTIONS TO READ RASTER ---------------------------------------------------------------------------------------------

def array_creation(data_array, value):
    """This function is to create an array of a coordinate

    -----------
    References:
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

    # Create zero array
    new_array = np.zeros((arr_x.shape[0], arr_y.shape[0]))

    # Get full number of x or y values
    if value == "x":
        for i in range(arr_x.shape[0]):
            for j in range(arr_y.shape[0]):
                new_array[i, j] = arr_x[i]

    else:
        for j in range(arr_y.shape[0]):
            for i in range(arr_x.shape[0]):
                new_array[i, j] = arr_y[j]


    return new_array


def xyz_array(data_array_func):
    """This function is to create an array from flowdepth raster including x, y, z values

    -----------
    References:
                None.
    -----------

    -----------
    Arguments:
                data_array_func:
                (array)
                                            Array of raster
    -----------

    -----------
    Returns:
                full_dataset:
                (array)
                                            A dataset contains 1D x, y, z array including padding
    -----------
    """

    # Create full number of values of x, y, z coordinates
    array_x = array_creation(data_array_func, 'x')
    array_y = array_creation(data_array_func, 'y')
    array_z = data_array_func.Band1.values

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

# END FUNCTION TO READ RASTER ------------------------------------------------------------------------------------------



# FUNCTION TO CREATE NEW INTERPOLATED LIDAR ----------------------------------------------------------------------------

def new_interpolation_data(data_array_func):
    """This function is to create new interpolated data by replace all current z values with new z values
    from 2D interpolation

    -----------
    References:
                None.
    -----------

    -----------
    Arguments:
                data_array_func:
                (array)
                                            Array of raster
    -----------

    -----------
    Returns:
                full_dataset:
                (array)
                                            A dataset contains 1D x, y, z array including padding
    -----------
    """

    # Make an array copy
    copy_full_data = data_array_func.copy()

    # Loop to replace current z values with new z values
    for i in range(data_array_func.shape[0]):
        copy_full_data[i, 2] = f(data_array_func[i, 0], data_array_func[i, 1])

    return copy_full_data

# END FUNCTION TO CREATE NEW INTERPOLATED LIDAR ------------------------------------------------------------------------


########################################################################################################################

# DEVELOP RIVER RASTER 1M ----------------------------------------------------------------------------------------------

crop_lidar(
    r"S:\\Lidar_example\\waikanae_part8.shp",
    r"S:\Lidar_example\points.laz",
    r"S:\Lidar_example\cropped_river_waikanae.laz",
    inside=True
)


# Interpolate at 1m resolution
instruction_interpolation = [
    r"S:\Lidar_example\cropped_river_waikanae.laz",
    {
        "type": "writers.gdal",
        "filename": r"S:\Lidar_example\river_raster.nc",
        "nodata": 0,
        "data_type": "float64",
        "output_type": "mean",
        "resolution": 1
    }
]

pipeline_interpolation = pdal.Pipeline(json.dumps(instruction_interpolation))
pipeline_interpolation.execute()

# Read raster and fill na
raster_river = rioxarray.open_rasterio(r"S:\Lidar_example\river_raster.nc")
new_raster_river = raster_river.rio.interpolate_na()

# Write filled raster into new raster
new_raster_river.rio.to_raster(r"S:\Lidar_example\waikanae_river_filled.nc")

# END DEVELOP RIVER RASTER 1M ------------------------------------------------------------------------------------------



# CROP ONLY RIVER AND CONVERT INTO LIDAR -------------------------------------------------------------------------------

# Read raster of nodata-filled full waikanae river
filled_raster_river = xarray.open_dataset(r"S:\Lidar_example\waikanae_river_filled.nc")

# Convert into 3D array (x, y, z)
full_data = xyz_array(filled_raster_river)

# Convert into LiDAR
writing_las_file(full_data)

# Crop only river
crop_lidar(
    r"S:\\Lidar_example\\waikanae_part8.shp",
    r"S:\Lidar_example\lidar_river_1m.laz",
    r"S:\Lidar_example\river_1m_crop.laz",
    inside=True
)

# Read the 1m resolution LiDAR that is just cropped
point_1m_class, point_1m_coord = read_las_file(r"S:\Lidar_example\river_1m_crop.laz")

# END CROP ONLY RIVER AND CONVERT INTO LIDAR ---------------------------------------------------------------------------





# DEVELOP FUNCTION FOR 2D INTERPOLATION WITH CUBIC --------------------------------------------------------------------

# Crop some parts of LiDAR for developing interpolation function
crop_lidar(
    fr"S:\Lidar_example\crop_before_1.shp",
    fr"S:\Lidar_example\cropped_river_waikanae.laz",
    fr"S:\Lidar_example\crop_before.laz",
    inside=True
)

# Read laz file
point_class, point_coord = read_las_file(r"S:\Lidar_example\crop_before.laz")

# Interpolation
f = interpolate.interp2d(point_coord[:, 0], point_coord[:, 1], point_coord[:, 2], kind="cubic")


# Crop 1m LiDAR
crop_lidar(
    fr"S:\Lidar_example\crop_before_1.shp",
    fr"S:\Lidar_example\river_1m_crop.laz",
    fr"S:\Lidar_example\crop_before_1m.laz",
    inside=True
)

point_1m_before_class, point_1m_before_coord = read_las_file(r"S:\Lidar_example\crop_before_1m.laz")

interpolated_before = new_interpolation_data(point_1m_before_coord)

# END DEVELOP FUNCTION FOR 2D INTERPOLATION WITH CUBIC -----------------------------------------------------------------




























# DEVELOP FUNCTION FOR 2D INTERPOLATION WITH CUBIC --------------------------------------------------------------------

# Crop some parts of LiDAR for developing interpolation function
num = 1
crop_lidar(
    fr"S:\Lidar_example\crop_{num}.shp",
    fr"S:\Lidar_example\crop_river_{num-1}.laz",
    fr"S:\Lidar_example\crop_river_{num}.laz",
    inside=False
)

# Read laz file
point_class, point_coord = read_las_file(r"S:\Lidar_example\crop_river_25.laz")

# Interpolation
f = interpolate.interp2d(point_coord[:, 0], point_coord[:, 1], point_coord[:, 2], kind="cubic")

# END DEVELOP FUNCTION FOR 2D INTERPOLATION WITH CUBIC -----------------------------------------------------------------


interpolated_data = new_interpolation_data(point_1m_coord)

writing_las_file(interpolated_data)


