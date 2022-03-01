# Prepare packages -----------------------------------------------------------------------------------------------------
import numpy as np                      # For data array manipulation
import pandas as pd                     # For reading data from text or csv files
import matplotlib.pyplot as plt         # For plotting for checking
# ----------------------------------------------------------------------------------------------------------------------


def discharge_data(discharge_file, time_file, resolution_func):
    """This function is to create discharge data following to resolution

    -----------
    References: http://www.bluesquarething.co.uk/geography/hydrograph.htm
                https://stackoverflow.com/questions/31789160/convert-select-columns-in-pandas-dataframe-to-numpy-array
    -----------

    -----------
    Arguments:
                discharge_file:
                (string)
                                        Path to the discharge file

                time_file:
                (string)
                                        Path to the time file

                resolution_func:
                (int or float)
                                        Resolution value in meter
    -----------

    -----------
    Returns:
               discharge_array:
               (array)
                                        A 2D array for discharge and time
    -----------

    """
    # Create 1D array of discharge
    reading_discharge = pd.read_csv(discharge_file, header=None)
    discharge = reading_discharge[0].to_numpy(dtype='float64')

    # Change discharge following to resolution
    discharge_res = discharge / resolution_func

    # Create 1D array of time
    reading_time = pd.read_csv(time_file, header=None)
    time = reading_time[0].to_numpy(dtype='float64')

    # Combine discharge_res and time into 2D array
    discharge_arr = np.column_stack((discharge_res, time))

    return discharge_arr


# Check data
if __name__ == '__main__':
    # Get discharge data
    discharge_array = discharge_data(10)

    # Plot hydrograph
    fig, ax = plt.subplots(figsize=(20, 20))
    ax.plot(discharge_array[:, 1], discharge_array[:, 0], "-o")
    plt.show()




