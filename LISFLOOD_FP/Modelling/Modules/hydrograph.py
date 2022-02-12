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
    # Small event
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

    # # Large event
    # discharge = np.array([
    #     0,
    #     20,
    #     100,
    #     200,
    #     400,
    #     100,
    #     200,
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
    #     1700,
    #     2000,
    #     2300,
    #     1700,
    #     1900,
    #     1550,
    #     1250,
    #     1000,
    #     900,
    #     700,
    #     400,
    #     200,
    #     100,
    #     50,
    #     100,
    #     20,
    #     10,
    #     5,
    #     0
    # ]) * 2

    # # Test case 1
    # discharge = np.array([
    #     0,
    #     1,
    #     3,
    #     4,
    #     10,
    #     25,
    #     50,
    #     40,
    #     20,
    #     10,
    #     5,
    #     2,
    #     0
    # ])

    # # For test case 2
    # discharge = np.array([
    #     0,
    #     2,
    #     5,
    #     8,
    #     9,
    #     10,
    #     12,
    #     15,
    #     10,
    #     8,
    #     5,
    #     2,
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

    # # For test case
    # time = np.array([
    #     0,
    #     10,
    #     20,
    #     30,
    #     40,
    #     50,
    #     60,
    #     70,
    #     80,
    #     90,
    #     100,
    #     110,
    #     120
    # ])




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




