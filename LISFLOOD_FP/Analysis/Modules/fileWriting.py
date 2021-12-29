# Preparing packages ---------------------------------------------------------------------------------------------------
from folder import *                            # For paths of sub-folders

# For handling raster
import numpy as np                              # For all calculation and data array/matrices manipulation
import rioxarray                                # For manipulating pixel values, spatial attributes, and raster files
import xarray                                   # For writing arrays into raster files
import rasterio                                 # For manipulating raster files (mainly crs and transformation)

# For handling polygons
from shapely.geometry import shape              # For converting geometry data into shapely geometry format
import geopandas as gpd                         # For manipulating geospatial data under polygon format
import pandas as pd                             # For manipulating dataframe
# ----------------------------------------------------------------------------------------------------------------------


# CSV WRITING ##########################################################################################################
def csv_generation(transformation_selection, dataset_func):
    """This function is to convert pandas dataframe into csv file to save the memory

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
                dataset_func:
                (pandas dataframe)
                                                 A dataset for being written into csv file
    -----------

    -----------
    Returns:
                None.
    -----------

    """
    # Set up the path for transformation_selection
    if transformation_selection == 'r':
        transformed = "rotated"
        csv_path = csv_rotation
    elif transformation_selection == 't':
        transformed = "translated"
        csv_path = csv_translation
    else:
        transformed = "combined"
        csv_path = csv_combination

    # Write into csv
    dataset_func.to_csv(fr"{csv_path}\\un_{transformed}_file_1.csv", index=False)
# END CSV WRITING ######################################################################################################


# RASTER WRITING #######################################################################################################
def raster_conversion(x_func, y_func, z_func):
    """This function is used to convert dimensions of array into raster array

    -----------
    References: https://stackoverflow.com/questions/41897544/make-a-contour-plot-by-using-three-1d-arrays-in-python
    -----------

    -----------
    Arguments:
                x_func:
                (array)
                                    1D x datasets
                y_func:
                (array)
                                    1D y datasets

                z_func:
                (array)
                                    1D z dataset
    -----------

    -----------
    Return:
                x_values_func:
                (array)
                                    Values of x under raster array (gridded array)
                y_values_func:
                (array)
                                    Values of y under raster array (gridded array)
                z_values_func:
                (array)
                                    Values of z under raster array (gridded array)
    -----------

    """
    # Gather x, y, z datasets into a pandas dataframe
    pd_dataframe = pd.DataFrame(dict(x=x_func, y=y_func, z=z_func))

    # Assign dataframe column names into variables
    xcol, ycol, zcol = 'x', 'y', 'z'

    # Sort dataframe according to x then y values
    pd_dataframe_sorted = pd_dataframe.sort_values(by=[xcol, ycol])

    # Getting values of x, y, z under raster array format
    # unique() function is used to remove duplicates
    x_values_func = pd_dataframe_sorted[xcol].unique()
    y_values_func = pd_dataframe_sorted[ycol].unique()
    z_values_func = pd_dataframe_sorted[zcol].values.reshape(len(x_values_func), len(y_values_func)).T

    return x_values_func, y_values_func, z_values_func


def raster_generation(transformation_selection, x_func, y_func, z_func, filename):
    """This function is used to write gridded array (x, y, z) into raster file

    -----------
    References: http://xarray.pydata.org/en/stable/user-guide/io.html
    -----------

    -----------
    Arguments:
                transformation_selection:
                (string)
                                            "r" means rotation
                                            "t" means translation
                                            "c" means combination
                x_func:
                (array)
                                            1D x datasets
                y_func:
                (array)
                                            1D y datasets

                z_func:
                (array)
                                            1D z dataset
                filename:
                (string)
                                            Name of raster file
    -----------

    -----------
    Returns:
                None.
    -----------

    """
    # Set up the path for transformation_selection
    if transformation_selection == 'r':
        transformed = "rotated"
        nc_raster_path_func = rotated_nc_raster_path
        raster_transformation_path = raster_rotation
    elif transformation_selection == 't':
        transformed = "translated"
        nc_raster_path_func = translated_nc_raster_path
        raster_transformation_path = raster_translation
    else:
        transformed = "combined"
        nc_raster_path_func = combined_nc_raster_path
        raster_transformation_path = raster_combination

    # Read original DEM raster without padding
    raster_origin_func = rioxarray.open_rasterio(fr"{nc_raster_path_func}\\generated_dem_no_padding.nc")

    # Round x and y to get duplicates
    x_round = np.round(x_func, 3)
    y_round = np.round(y_func, 3)

    # Get x, y, z values under raster array (gridded array)
    x_values_func, y_values_func, z_values_func = raster_conversion(x_round, y_round, z_func)

    # Write x, y, z values into raster array
    raster_array = xarray.DataArray(
        data=z_values_func,
        dims=['y', 'x'],
        coords={
            'x': (['x'], x_values_func),
            'y': (['y'], y_values_func)
        },
        attrs=raster_origin_func.attrs
    )

    # Set up crs and nodata
    raster_array.rio.write_crs("epsg:2193", inplace=True)
    raster_array.rio.write_nodata(-999, inplace=True)

    # Write into raster file (tiff)
    raster_array.rio.to_raster(fr"{raster_transformation_path}\\{filename}_un{transformed}_flowdepth_raster.nc")
# END RASTER WRITING ###################################################################################################



# POLYGON WRITING ######################################################################################################
def polygon_conversion(transformation_selection, dataset_func, column):
    """This function is to convert dataset into polygon (to write into shapefile)

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
                dataset_func:
                (pandas dataframe)
                                            A dataset that needs converting to polygon
                column:
                (string)
                                            Column name that needs converting
    -----------

    -----------
    Returns:
                poly_dataframe:
                (pandas dataframe)
                                            A geopandas dataframe
    -----------

    """
    # Set up the path for transformation_selection
    if transformation_selection == 'r':
        nc_raster_path_func = rotated_nc_raster_path
    elif transformation_selection == 't':
        nc_raster_path_func = translated_nc_raster_path
    else:
        nc_raster_path_func = combined_nc_raster_path

    # Convert 0 values into NaN
    dataset_func.loc[dataset_func[f"{column}"] == -999, [f"{column}"]] = np.nan

    # Read original DEM raster without padding
    poly_origin = rasterio.open(fr"{nc_raster_path_func}\\generated_dem_no_padding.nc")

    # Get information from original DEM raster
    poly_origin_array = poly_origin.read(1)
    poly_origin_transform = poly_origin.transform
    poly_origin_crs = poly_origin.crs

    # Extract parameters: id and depth
    id_pixels = np.arange(poly_origin_array.size).reshape(poly_origin_array.shape)
    value_list = dataset_func[f"{column}"].tolist()

    # Vectorise features
    vectors = rasterio.features.shapes(source=id_pixels.astype(np.int16), transform=poly_origin_transform)

    vectors_list = list(vectors)

    # Get geometry
    poly_geometry = [shape(polygon) for polygon, value in vectors_list]

    # Get id
    id_poly = [id_val for id_poly, id_val in vectors_list]

    # Create database
    poly_database = {'id': id_poly,
                     f"{column}": value_list}
    poly_dataframe = gpd.GeoDataFrame(data=poly_database,
                                      geometry=poly_geometry,
                                      crs=poly_origin_crs)

    return poly_dataframe


def polygon_generation(transformation_selection, dataset_func, column):
    """This function is to write dataframe into shape file

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
                dataset_func:
                (pandas dataframe)
                                                A dataset that needs converting to polygon
                column:
                (string)
                                                Column name that needs converting
    -----------

    -----------
    Returns:
                None.
    -----------

    """
    # Set up the path for transformation_selection
    if transformation_selection == 'r':
        transformed = "rotated"
        polygon_transformation_path = polygon_rotation
    elif transformation_selection == 't':
        transformed = "translated"
        polygon_transformation_path = polygon_translation
    else:
        transformed = "combined"
        polygon_transformation_path = polygon_combination

    # Get polygon dataframe
    polygon_dataframe = polygon_conversion(transformation_selection, dataset_func, column)

    # Write into shape file
    polygon_dataframe.to_file(fr"{polygon_transformation_path}\\{column}_un{transformed}_flowdepth_polygon.shp",
                              driver='ESRI Shapefile')
# END POLYGON WRITING ##################################################################################################