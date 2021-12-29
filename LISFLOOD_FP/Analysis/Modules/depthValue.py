# Prepare packages -----------------------------------------------------------------------------------------------------
from folder import *                                    # For paths of sub-folders
import time                                             # For timing steps

# For handling polygons
import fiona                                            # For reading shape file under GeoJSON format
from shapely.geometry import shape                      # For converting geometry data into shapely geometry format
from rtree import index                                 # For using spatial index function
# ----------------------------------------------------------------------------------------------------------------------


def get_flowdepth_value(transformation_selection, number_simulation, time_extract_func):
    """This function is to extract the flowdepth at given coordinates from any shape file

    -----------
    References: https://gis.stackexchange.com/questions/102933/more-efficient-spatial-join-in-python-without-qgis-arcgis-postgis-etc/103066#103066
                https://gis.stackexchange.com/questions/121469/get-shapefile-polygon-attribute-value-at-a-specific-point-using-python-e-g-via
                https://stackoverflow.com/questions/59030022/checking-whether-point-is-within-polygon-returns-wrong-results-in-shapely
                https://gis.stackexchange.com/questions/119919/maximizing-code-performance-for-shapely
                https://gis.stackexchange.com/questions/42931/rtree-python-polygon-index
                https://gis.stackexchange.com/questions/227474/rtree-spatial-index-does-not-result-in-faster-intersection-computation

                https://rtree.readthedocs.io/en/latest/tutorial.html

                https://sgillies.net/2014/01/18/getting-shapes-of-raster-features-with-rasterio.html
                https://stackoverflow.com/questions/20474549/extract-points-coordinates-from-a-polygon-in-shapely
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
                flowdepth:
                (list)
                                                A list of flowdepth values
    -----------

    """

    # Set up the path for transformation_selection
    if transformation_selection == 'r':
        transformed = "rotated"
        untransformed_path = unrotated_path
    elif transformation_selection == 't':
        transformed = "translated"
        untransformed_path = untranslated_path
    else:
        transformed = "combined"
        untransformed_path = uncombined_path

    # A list stores flowdepth values
    flowdepth_list = []

    # Input polygons and points files
    polygons_input = fiona.open(
        fr"{untransformed_path}\shapefile_{number_simulation}\flowdepth_un{transformed}_{number_simulation}_at_{time_extract_func}.geojson")
    points_input = fiona.open(f"{untransformed_path}\\shapefile_Point\\Point_clip.shp")

    # Read polygons and points files into geometry lists
    values_polygons_input = [values_poly for values_poly in polygons_input]
    values_points_input = [values_point for values_point in points_input]

    # Construct r-tree spatial index
    indexes = index.Index()

    # Insert boundary records obtained from polygons geometry using shapely shape() function
    # In other words, create boundaries based on the coordinates gained from
    # polygon geometry under shapely format (left, bottom, right, top)
    for id_poly, polygons in enumerate(values_polygons_input):
        indexes.insert(id_poly, shape(polygons['geometry']).bounds)

    # Iterate through each point and polygon
    # to find out the spatial index of polygon that contains the point.
    # In other words, checking if each point fits in the given boundaries
    # (boundaries here are just a gather of rectangle boundaries created above)
    # After that, check if the polygon with that spatial index contains the point.
    # The flowdepth is extracted from that polygon
    # (assumed that all coordinates within that polygon have the same flowdepth)

    # Loop 1: Convert point under GeoJSON into shapely geometry
    for id_point, points in enumerate(values_points_input):
        each_point = shape(points['geometry'])

        # Loop 2: Check if the point is in any boundaries and polygons
        for id_intersection in indexes.intersection(each_point.coords[0]):
            # The above line will return the list of spatial indexes that have intersection
            # between the point and any boundaries.

            # Check if the point is in the polygon identified by the spatial index.
            # within() is used to find out the polygon the point lies in
            if each_point.within(shape(values_polygons_input[id_intersection]['geometry'])):
                flowdepth_value = values_polygons_input[id_intersection]['properties']['depth']
                flowdepth_list.append(flowdepth_value)

    return flowdepth_list


def run_depth_value_extraction(set_of_simulation, clipped_dataset_func,
                               transformation_selection, time_extract_func):
    """This function is to get the list of all transformed values from all simulations

    -----------
    References:
                None.
    -----------

    -----------
    Arguments:
                set_of_simulation:
                (list)
                                                An 3D array of simulations (angle, x, y)
                clipped_dataset_func:
                (dictionary)
                                                A dictionary already contains x, y coordinate values
                transformation_selection:
                (string)
                                                "r" means rotation
                                                "t" means translation
                                                "c" means combination
                time_extract_func:
                (int)
                                                Amount of time that flood model predicted
    -----------

    -----------
    Returns:
                None.
    -----------

    """
    # Run simulations to collect depth values
    for order_number in set_of_simulation:
        # Get values
        angle_val = order_number[0]
        x_val = order_number[1]
        y_val = order_number[2]

        # Set the name of each transformation
        combination_number = f"angle_{angle_val}_x_{x_val}_y_{y_val}"

        # Print the header of the first running
        print("{0} angle: {1} and x: {2} and y: {3} {0}".format("-" * 30, angle_val, x_val, y_val))

        # Start timing the depth values extraction process -------------------------------------------------------------
        start_depth = time.time()
        # Append depth values into the dictionary
        clipped_dataset_func[f"{combination_number}"] = get_flowdepth_value(transformation_selection,
                                                                            combination_number, time_extract_func)
        # End timing the depth values extraction process ---------------------------------------------------------------
        end_depth = time.time()

        # PRINT --------------------------------------------------------------------------------------------------------
        # Get the running time
        depth_time = end_depth - start_depth

        # Convert second to minute
        depth_minute = depth_time / 60

        # Get the name
        depth_name = "* Running time of depth value extraction:"

        # Print
        print('{0:<60}{1:>10}'.format(depth_name, depth_minute))

        # Print separator
        print("-" * 90)
        print("\n")