# Preparing packages ---------------------------------------------------------------------------------------------------
from folder import *                            # For paths of sub-folders

# For handling raster
import numpy as np                              # For all calculation and data array/matrices manipulation
import rioxarray                                # For manipulating pixel values, spatial attributes, and raster files
import xarray                                   # For writing arrays into raster files
import rasterio                                 # For manipulating raster files (mainly crs and transformation)
import rasterio.features                        # For manipulating raster features
from datasetClipping import array_creation      # For getting data array shape from raster data


# For handling polygons
from shapely.geometry import shape              # For converting geometry data into shapely geometry format
import geopandas as gpd                         # For manipulating geospatial data under polygon format
import pandas as pd                             # For manipulating dataframe
from pyogrio import write_dataframe             # For writing out shape file (twice faster than geopandas)

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
def raster_conversion(transformation_selection, x_func, y_func, z_func, rectangle=True, switch=False):
    """This function is used to convert dimensions of array into raster array

    -----------
    References: https://stackoverflow.com/questions/41897544/make-a-contour-plot-by-using-three-1d-arrays-in-python
                https://stackoverflow.com/questions/41815079/pandas-merge-join-two-data-frames-on-multiple-columns
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
                rectangle:
                (boolean)
                                            To specify that if the data shape is rectangle or not.
                                            True is rectangle (default)
                                            Fasle is not rectangle
                switch:
                (boolean)
                                            To switch the row and column of array
                                            True is switching
                                            False is not switching (default)
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
    # Set up the path for transformation_selection
    if transformation_selection == 'r':
        nc_raster_path_func = rotated_nc_raster_path
    elif transformation_selection == 't':
        nc_raster_path_func = translated_nc_raster_path
    else:
        nc_raster_path_func = combined_nc_raster_path

    # Import a standard raster to extract its attributes
    raster_origin_func = rioxarray.open_rasterio(fr"{nc_raster_path_func}\\generated_dem_no_padding.nc")

    # Gather x, y, z datasets into a pandas dataframe (DATAFRAME OF ACTUAL DATA)
    pd_dataframe_actual = pd.DataFrame(dict(x=x_func, y=y_func, z=z_func))

    # Assign dataframe column names into variables
    xcol, ycol, zcol = 'x', 'y', 'z'

    # Sort dataframe according to x then y values
    pd_dataframe_sorted = pd_dataframe_actual.sort_values(by=[xcol, ycol])

    # If the data shape (array shape of flatten x, y, z) is rectangle, we can just convert data into x, y, z format of
    # raster generation. If the data shape is not rectangle, we should create a background (or called fake data) with a
    # specific value (normally 0 will be chosen), then merge it with in-rectangular data. In specific, we should create
    # two datasets, one dataset is the array of x, y, z representing for the shape of the data, and the other dataset
    # is the array of x, y, z representing for the shape of background (should be rectangular). We then convert these
    # datasets under numpy array format into pandas dataframes with three columns x, y, z. We then merge them based
    # on x and y columns into a new dataframe where no-data areas in in-rectangular data are filled with specific
    # values (0 maybe)

    if rectangle:
        # Getting values of x, y, z under raster array format
        # unique() function is used to remove duplicates
        x_values_func = pd_dataframe_sorted[xcol].unique()
        y_values_func = pd_dataframe_sorted[ycol].unique()
        z_values_func = pd_dataframe_sorted[zcol].values.reshape(len(x_values_func), len(y_values_func)).T

    else:
        # CREATE DATAFRAME OF BACKGROUND DATA
        # Getting values of x, y, z under raster array format for background data
        # unique() function is used to remove duplicates
        x_background = pd_dataframe_sorted[xcol].unique()
        y_background = pd_dataframe_sorted[ycol].unique()
        z_background = np.full((len(x_background), len(y_background)), -999).T

        # Create a raster data from above background data
        raster_background = xarray.DataArray(
            data=z_background,
            dims=['y', 'x'],
            coords={
                'x': (['x'], x_background),
                'y': (['y'], y_background)
            },
            attrs=raster_origin_func.attrs
        )

        # Getting x, y ,z array from raster background
        array_x_background = array_creation(raster_background, 'x', switch)
        array_y_background = array_creation(raster_background, 'y', switch)
        array_z_background = raster_background.values

        # Flatten all x, y, z
        flatten_x_background = array_x_background.flatten()
        flatten_y_background = array_y_background.flatten()
        flatten_z_background = array_z_background.flatten()

        # Create 3d array based on flattened x, y, z
        full_dataset_background = np.vstack((flatten_x_background,
                                             flatten_y_background,
                                             flatten_z_background)).transpose()

        # Dataframe of background data
        pd_dataframe_background = pd.DataFrame(dict(x=full_dataset_background[:, 0],
                                                y=full_dataset_background[:, 1],
                                                z=full_dataset_background[:, 2]))

        # MERGING BACKGROUND DATA & ACTUAL DATA
        merged_data = pd.merge(pd_dataframe_background, pd_dataframe_actual,
                               how='left', left_on=['x', 'y'], right_on=['x', 'y'])

        # COPY THE MERGED DATA & FILL NAN VALUES WITH A SPECIFIC VALUE (0 maybe)
        copy_merged_data = merged_data[['x', 'y', 'z_y']].copy()
        copy_merged_data['z_y'] = copy_merged_data['z_y'].fillna(-999)

        # SORT DATA BEFORE SEPARATE THEM INTO X, Y, Z
        # Name the columns
        xcol_now, ycol_now, zcol_now = 'x', 'y', 'z_y'

        # Sort them
        copy_merged_data_sorted = copy_merged_data.sort_values(by=[xcol_now, ycol_now])

        # WRITE DATA INTO RASTER FORMAT
        x_values_func = copy_merged_data_sorted['x'].unique()
        y_values_func = copy_merged_data_sorted['y'].unique()
        z_values_func = copy_merged_data_sorted['z_y'].values.reshape(len(x_values_func), len(y_values_func)).T

    return x_values_func, y_values_func, z_values_func


def raster_generation(transformation_selection, x_func, y_func, z_func, filename, save_path=None, rectangle=True,
                      switch=False):
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
                save_path:
                (string)
                                            If None, save rasters of statistical results.
                                            If specified, save rasters to specified path.
                rectangle:
                (boolean)
                                            To specify that if the data shape is rectangle or not.
                                            True is rectangle (default)
                                            False is not rectangle
                switch:
                (boolean)
                                            To switch the row and column of array
                                            True is switching
                                            False is not switching (default)
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
    x_values_func, y_values_func, z_values_func = raster_conversion(transformation_selection,
                                                                    x_round, y_round, z_func,
                                                                    rectangle, switch)

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

    # Write into raster file (netcdf)
    if save_path is None:
        raster_array.rio.to_raster(fr"{raster_transformation_path}\\{filename}_un{transformed}_flowdepth_raster.nc")

    else:
        raster_array.rio.to_raster(fr"{save_path}\\{filename}.nc")
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


def polygon_generation(transformation_selection, dataset_func, column, save_path):
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
                save_path:
                (string)
                                                If None, save rasters of statistical resutls.
                                                If specify, save rasters to specified path.
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
    if save_path is None:
        write_dataframe(polygon_dataframe,
                        fr"{polygon_transformation_path}\\{column}_un{transformed}_flowdepth_polygon.shp",
                        driver='GeoJSON')
    else:
        write_dataframe(polygon_dataframe,
                        fr"{save_path}\\{column}_un{transformed}_flowdepth_polygon.shp",
                        driver='GeoJSON')

# END POLYGON WRITING ##################################################################################################