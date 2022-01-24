# Prepare packages -----------------------------------------------------------------------------------------------------
import numpy as np
import matplotlib.pyplot as plt
# ----------------------------------------------------------------------------------------------------------------------


def discharge_data(resolution_func):
    """This function is to create discharge data following to resolution

    -----------
    References: http://www.bluesquarething.co.uk/geography/hydrograph.htm
    -----------

    -----------
    Arguments:
               resolution_func:
               (int or float)
                                            resolution value in meter
    -----------

    -----------
    Returns:
               discharge_array:
               (array)
                                    A 2D array for discharge and time
    -----------

    """
    # Create 1D array of discharge
    discharge = np.array([
        5,
        5,
        10,
        10,
        20,
        20,
        60,
        100,
        200,
        250,
        190,
        90,
        40,
        20,
        50,
        120,
        200,
        500,
        520,
        550,
        900,
        980,
        1000,
        1300,
        1000,
        950,
        800,
        550,
        220,
        100,
        50,
        20,
        10,
        5,
        2,
        1,
        0
    ])

    # discharge = np.array([
    #     0,
    #     50,
    #     300,
    #     400,
    #     500,
    #     700,
    #     450,
    #     400,
    #     350,
    #     400,
    #     650,
    #     750,
    #     850,
    #     900,
    #     1100,
    #     1230,
    #     1450,
    #     1750,
    #     2200,
    #     2350,
    #     2050,
    #     1850,
    #     1900,
    #     1550,
    #     1250,
    #     1230,
    #     1150,
    #     1150,
    #     700,
    #     550,
    #     350,
    #     300,
    #     200,
    #     50,
    #     30,
    #     10,
    #     0
    # ])


    # Change discharge following to resolution
    discharge_res = discharge / resolution_func

    # Create 1D array of time
    time = np.array([
        0,
        200,
        400,
        600,
        800,
        1000,
        1200,
        1400,
        1600,
        1800,
        2000,
        2200,
        2400,
        2600,
        2800,
        3000,
        3200,
        3400,
        3600,
        3800,
        4000,
        4200,
        4400,
        4600,
        4800,
        5000,
        5200,
        5400,
        5600,
        5800,
        6000,
        6200,
        6400,
        6600,
        6800,
        7000,
        7200
    ])

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




