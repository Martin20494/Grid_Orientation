# Prepare packages -----------------------------------------------------------------------------------------------------
from folder import *                                # For paths of sub-folders
import glob                                         # For getting all files' names in a certain path
import rioxarray                                    # For manipulating spatial data under xarray array format
# ----------------------------------------------------------------------------------------------------------------------


# Create function to extract the flowdepth from BG_Flood output
def flowdepth_extraction(transformation_selection, number_simulation, time_extract_func):
    """This function is to extract the flowdepth from LISFLOOD-FP output

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
                number_simulation:
                (int)
                                            Ordinal number of simulation

                time_extract_func:
                (int)
                                            Amount of time that flood model predicted
    -----------

    -----------
    Returns:
                None.
    -----------

    """
    # Set up the path for transformation_selection
    if transformation_selection == 'r':
        output_path = rotated_FPoutput_path
        flowdepth_path = rotated_flowdepth
        transformed = "rotated"
    elif transformation_selection == 't':
        output_path = translated_FPoutput_path
        flowdepth_path = translated_flowdepth
        transformed = "translated"
    else:
        output_path = combined_FPoutput_path
        flowdepth_path = combined_flowdepth
        transformed = "combined"

    # Extract required flowdepth file
    raster_asc = rioxarray.open_rasterio(glob.glob(fr"{output_path}\\{transformed}_{number_simulation}\\*{time_extract_func}.wd")[0])

    # Add crs
    new_raster = raster_asc.rio.write_crs(2193)

    # Write to flowdepth folder
    new_raster.rio.to_raster(fr"{flowdepth_path}\\flowdepth_{transformed}_{number_simulation}_at_{time_extract_func}.nc")
