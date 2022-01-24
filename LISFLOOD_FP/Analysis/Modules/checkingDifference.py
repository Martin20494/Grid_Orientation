# Prepare packages -----------------------------------------------------------------------------------------------------
import seaborn as sns                                           # For plotting density and frequency
import numpy as np                                              # For data handling
# ----------------------------------------------------------------------------------------------------------------------

def get_diff(full_data_func):
    """This function is to get the differences in the outputs of 0 and 90 degrees simulations

    -----------
    References:
                None.
    -----------

    -----------
    Arguments:
                full_data_func:
                (pandas dataframe)
                                        A dataframe contains simulations of 0 and 90 degrees
    -----------

    -----------
    Returns:
                diff_list_func:
                (list)
                                        A list of difference in values of each cell between simulations of 0 and 90
                                        degrees' depth values
    -----------

    """
    # Set up the list of difference
    diff_list_func = []

    # Collect different values
    for each_cell in range(full_data_func.shape[0]):
        # Calculate the difference
        diff_value = full_data_func['angle_90_x_0_y_0'][each_cell] - full_data_func['angle_0_x_0_y_0'][each_cell]

        # Check if the difference is zero means if simulations of 0 and 90 have different values to put into the list
        if diff_value != 0:
            diff_list_func.append(diff_value)

    return diff_list_func



def plot_diff(diff_list_func, axis_func, resolution_func, interval_func):
    """This function is to plot the histogram of the differences between the simulations of 0 and 90 degrees' depth
    values

    -----------
    References:
                https://stackoverflow.com/questions/12608788/changing-the-tick-frequency-on-x-or-y-axis-in-matplotlib
    -----------

    -----------
    Arguments:
                diff_list_func:
                (list)
                                        A list of difference in values of each cell between simulations of 0 and 90
                                        degrees' depth values
                axis_func:
                (axis in matplotlib)
                                        The order of axis in matplotlib figure
                resolution_func:
                (int or float)
                                        Resolution value in meter
                interval_func:
                (int)
                                        Interval value for y axis
    -----------

    -----------
    Returns:
                None.
    -----------

    """
    # Get the minimum and maximum of the differences
    print("Minimum different water depth value:", min(diff_list_func))
    print("Maximum different water depth value:", max(diff_list_func))

    # Draw histogram
    sns.distplot(
        diff_list_func,
        bins=100,
        kde=False,
        ax=axis_func,
        hist_kws={
            "histtype": "bar",
            "linewidth": 2,
            "alpha": 1,
            "color": "b",
            "edgecolor": "black"
        }
    )

    # Add titles and labels
    title = f'Different water depth values\nof each cell between 0 and 90 degree,\n'
    title += fr"resolution = {resolution_func} meters"
    axis_func.set_title(title, pad=20, fontsize=25, fontweight='bold')
    axis_func.set_xlabel("Different water depth values (m)", labelpad=38, fontsize=20)
    axis_func.set_ylabel("Frequency (number of cells)", labelpad=38, fontsize=20)

    # Increase figures in axes' labels
    for item in (axis_func.get_xticklabels() + axis_func.get_yticklabels()):  # For x, y ticks' labels
        item.set_fontsize(15)

    # Increase the size of ticks
    axis_func.tick_params(direction='out', length=8, pad=10)

    # Get integers for ylim
    # https://stackoverflow.com/questions/12608788/changing-the-tick-frequency-on-x-or-y-axis-in-matplotlib (for
    # changing interval in ylim
    start, end = axis_func.get_ylim()
    axis_func.yaxis.set_ticks(np.arange(start, end, interval_func))

    # Remove grid background lines (including x, y lines)
    axis_func.grid(False)
    axis_func.spines['top'].set_visible(False)
    axis_func.spines['right'].set_visible(False)
    axis_func.spines['bottom'].set_visible(False)
    axis_func.spines['left'].set_visible(False)


def difference_information(full_data_func, diff_list_func):
    """This function is to get information from different list

    -----------
    References:
                None.
    -----------

    -----------
    Arguments:
                full_data_func:
                (pandas dataframe)
                                        A dataframe contains simulations of 0 and 90 degrees
                diff_list_func:
                (list)
                                        A list of difference in values of each cell between simulations of 0 and 90
                                        degrees' depth values
    -----------

    -----------
    Returns:
                None.
    -----------

    """
    # Change diff_list_func into diff_array_func
    diff_array_func = np.array(diff_list_func).astype('float64')

    # Get titles
    title1 = "Total cells:"
    title2 = "Total cells having 0-90 different values:"
    title3 = "Total cells having absolute 0-90 different values > 0.001:"
    title4 = "Minimum:"
    title5 = "Maximum:"

    # Print results
    print('\033[1m{0:<60}\033[0m{1:<3}'.format(title1, full_data_func.shape[0]))          # For integers, using '<5' for left
                                                                            # alignment
    print('\033[1m{0:<60}\033[0m{1:<3} ({2:.4f} %)'.format(
        title2,
        diff_array_func.shape[0],
        diff_array_func.shape[0]/full_data_func.shape[0]*100
    ))                                                                      # For integers, using '<5' for left
                                                                            # alignment
    print('\033[1m{0:<60}\033[0m{1:<3} ({2:.4f} %)'.format(
        title3,
        diff_array_func[np.abs(diff_array_func) > 0.001].shape[0],
        diff_array_func[np.abs(diff_array_func) > 0.001].shape[0]/full_data_func.shape[0]*100
    ))                                                                      # For integers, using '<5' for using left
                                                                            # alignment
    print('\033[1m{0:<60}\033[0m{1:>3}'.format(title4, np.min(diff_array_func)))          # For float, using '>5' for left alignment
    print('\033[1m{0:<60}\033[0m{1:>3}'.format(title5, np.max(diff_array_func)))          # For float, using '>5' for left alignment