# Prepare packages -----------------------------------------------------------------------------------------------------
from folder import *                                             # For paths of sub-folders
# ----------------------------------------------------------------------------------------------------------------------


def save_plot(transformation_selection, axis_func, fig_func, plot_name, extension, x_portion, y_portion):
    """This function is to save the plots

    -----------
    References: https://stackoverflow.com/questions/4325733/save-a-subplot-in-matplotlib
    -----------

    -----------
    Arguments:
                transformation_selection:
                (string)
                                            "r" means rotation
                                            "t" means translation
                                            "c" means combination
                axis_func:
                (axis in matplotlib)
                                            The ordinal number of plot in subplot
                fig_func:
                (figure in matplotlib)
                                            Figure in subplot matplotlib
                plot_name:
                (string)
                                            Name of the plot
                extension:
                (string)
                                            "png" or "pdf"
                x_portion:
                (int)
                                            Portion of width of plot in the saved filed (%).
                                            For example, if 20%, input 2
                y_portion:
                (int)
                                            Portion of length of plot in the saved filed (%).
                                            For example, if 20%, input 2
    -----------

    -----------
    Returns:
                None.
    -----------

    """
    if transformation_selection == "r":
        saving_path = plot_rotation
    elif transformation_selection == "t":
        saving_path = plot_translation
    else:
        saving_path = plot_combination

    # Create a parameter to control the size of the plot when saving
    extent = axis_func.get_window_extent().transformed(fig_func.dpi_scale_trans.inverted())

    # Put the parameter to the plot
    fig_func.savefig(fr"{saving_path}\\{plot_name}.{extension}", bbox_inches=extent.expanded(float(f'1.{x_portion}'),
                                                                                             float(f'1.{y_portion}')))