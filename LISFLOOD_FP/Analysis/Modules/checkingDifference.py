# Prepare packages -----------------------------------------------------------------------------------------------------
import seaborn as sns                                           # For plotting density and frequency
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



def plot_diff(diff_list_func, axis_func):
    """This function is to plot the histogram of the differences between the simulations of 0 and 90 degrees' depth
    values

    -----------
    References:
                None.
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
    axis_func.set_title("Different water depth values of each cell between 0 and 90 degree", pad=20, fontsize=25,
                 fontweight='bold')
    axis_func.set_xlabel("Different water depth values (m)", labelpad=38, fontsize=20)
    axis_func.set_ylabel("Frequency (number of cells)", labelpad=38, fontsize=20)

    # Increase figures in axes' labels
    for item in (axis_func.get_xticklabels() + axis_func.get_yticklabels()):  # For x, y ticks' labels
        item.set_fontsize(15)

    # Increase the size of ticks
    axis_func.tick_params(direction='out', length=8, pad=10)

    # Remove grid background lines (including x, y lines)
    axis_func.grid(False)
    axis_func.spines['top'].set_visible(False)
    axis_func.spines['right'].set_visible(False)
    axis_func.spines['bottom'].set_visible(False)
    axis_func.spines['left'].set_visible(False)