# Prepare packages -----------------------------------------------------------------------------------------------------

# For controlling the paths
from folder import *                              # For paths of sub-folders

# For data handling generally
import numpy as np                                # For all calculation and data array/matrices manipulation

# For raster and vector handling
import rasterio                                   # For reading and manipulating spatial data
import rasterio.features                          # For vectorising features in array
from fileWriting import raster_generation         # For generating raster


# For polygon handling
import geopandas as gpd                           # For manipulating shape files
from shapely.geometry import shape                # For manipulating spatial information (geometry) under GeoJSON format
from shapely.ops import unary_union               # For combining all polygons into one
from pyogrio import write_dataframe               # For writing out shape file (twice faster than geopandas)
# ----------------------------------------------------------------------------------------------------------------------

def one_polygon_raster_generation(transformation_selection, dataset_func, filter_rate_func):
    """This function is to convert dataset into raster which will be used for converting into one polygon

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
                filter_rate_func:
                (float or int)
                                            The rate at which the depth values will be ignored
    -----------

    -----------
    Returns:
                None.
    -----------

    """
    # Set up the path for transformation_selection
    if transformation_selection == 'r':
        transformed = "rotated"
        one_poly_path = one_polygon_rotation
    elif transformation_selection == 't':
        transformed = "translated"
        one_poly_path = one_polygon_translation
    else:
        transformed = "combined"
        one_poly_path = one_polygon_combination

    # Create path to store file
    oneraster_path = f"{one_poly_path}\\oneRaster"
    pathlib.Path(oneraster_path).mkdir(parents=True, exist_ok=True)

    for num in range(len(dataset_func.columns) - 2):    # Minus 2 due to two columns of coordinates x and y
        # Convert to -999
        col = dataset_func.columns[num + 2]
        dataset_func.loc[dataset_func[col] < filter_rate_func, [col]] = -999
        
        # Generate raster
        raster_generation(
            transformation_selection,
            dataset_func['x_coord'],
            dataset_func['y_coord'],
            dataset_func[f'{dataset_func.columns[num + 2]}'],
            f"oneRaster_un{transformed}_{dataset_func.columns[num + 2]}",
            oneraster_path
        )



def one_polygon_conversion(transformation_selection, dataset_func, column, filter_rate_func):
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
                filter_rate_func:
                (float or int)
                                            The rate at which the depth values will be ignored
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
        transformed = "rotated"
        one_poly_path = one_polygon_rotation
    elif transformation_selection == 't':
        transformed = "translated"
        one_poly_path = one_polygon_translation
    else:
        transformed = "combined"
        one_poly_path = one_polygon_combination

    # Read original DEM raster without padding
    raster_poly = rasterio.open(fr"{one_poly_path}\\oneRaster\\oneRaster_un{transformed}_{column}.nc")

    # Get information from original DEM raster
    raster_array = raster_poly.read(1)
    raster_transform = raster_poly.transform
    raster_crs = raster_poly.crs

    # Extract parameters: id and depth
    id_pixels = np.arange(raster_array.size).reshape(raster_array.shape)
    value_list = dataset_func[f"{column}"].tolist()

    # Vectorise features
    vectors = rasterio.features.shapes(source=id_pixels.astype(np.int16), transform=raster_transform)

    vectors_list = list(vectors)

    # Get geometry
    poly_geometry = [shape(polygon) for polygon, value in vectors_list]

    # Get id
    id_poly = [id_val for id_poly, id_val in vectors_list]

    # Create database
    poly_database = {'id': id_poly,
                     "depth": value_list}
    poly_dataframe = gpd.GeoDataFrame(data=poly_database,
                                      geometry=poly_geometry,
                                      crs=raster_crs)

    # Select only FLOODED area by using the filter rate
    flood_polygons_list = poly_dataframe[poly_dataframe['depth'] >= filter_rate_func]['geometry'].tolist()

    # Combine all polygons into one using function unary_union of shapely package
    one_flood_polygon = gpd.GeoSeries(unary_union(flood_polygons_list))

    # Write that one polygon again into geopandas dataframe
    one_poly_gdf = gpd.GeoDataFrame(data={"id": [1]},
                                    geometry=one_flood_polygon.tolist(),
                                    crs=2193)

    # Create path
    new_onepoly_path = f"{one_poly_path}\\onePoly"
    pathlib.Path(new_onepoly_path).mkdir(parents=True, exist_ok=True)

    write_dataframe(one_poly_gdf,
                    fr"{new_onepoly_path}\\onepoly_un{transformed}_{column}.geojson",
                    driver="GeoJSON")


def one_polygon_generation(transformation_selection, dataset_func, filter_rate_func):
    """This function is to write dataframe into one-polygon shape file

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
                filter_rate_func:
                (float or int)
                                                The rate at which the depth values will be ignored
    -----------

    -----------
    Returns:
                None.
    -----------

    """
    for num in range(len(dataset_func.columns) - 2):
        one_polygon_conversion(
            transformation_selection,
            dataset_func,
            dataset_func.columns[num + 2],
            filter_rate_func
        )