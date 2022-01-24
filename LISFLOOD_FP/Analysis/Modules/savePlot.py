# Prepare packages -----------------------------------------------------------------------------------------------------
from folder import *                                             # For paths of sub-folders
# ----------------------------------------------------------------------------------------------------------------------


def save_plot(transformation_selection, fig_func, plot_name, extension, dpi_func, addition_path=""):
    """This function is to save the plots

    -----------
    References: https://stackoverflow.com/questions/4325733/save-a-subplot-in-matplotlib
                https://stackoverflow.com/questions/16032389/pad-inches-0-and-bbox-inches-tight-makes-the-plot-smaller-than-declared-figsiz
                https://stackoverflow.com/questions/10041627/how-to-make-savefig-save-image-for-maximized-window-instead-of-default-size
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
                dpi_func:
                (int or float)
                                            The resolution in dots per inch. Please visit here for more information
                                            https://matplotlib.org/stable/api/_as_gen/matplotlib.pyplot.savefig.html
                addition_path:
                (string)
                                            Path for saving plots
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
    elif transformation_selection == "c":
        saving_path = plot_combination
    else:
        saving_path = addition_path

    fig_func.savefig(fr"{saving_path}\\{plot_name}.{extension}", bbox_inches='tight', dpi=dpi_func)