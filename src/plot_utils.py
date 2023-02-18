""" Plotting utilities """ 

import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

def plot_thickness_colour(thickness, film_colour_sRGB, name="", title="", save=0):
    """ Plot showing the variation of colour with thickness """
    # Creates colourmap from RGB values
    custom_cmap = mpl.colors.ListedColormap(film_colour_sRGB)
    # sets min and max of colour bar labels
    norm = mpl.colors.Normalize(vmin=np.min(thickness), vmax=np.max(thickness))
    
    fs_labels = 12
    fs_ticks = 11
    width = 12.9 / 2.54  # inches
    fig, ax = plt.subplots(figsize=(width, 0.30 * width), constrained_layout=True,
                            dpi=123.6)
    bar = mpl.colorbar.ColorbarBase(ax, cmap = custom_cmap, norm = norm, 
                                    orientation = 'horizontal')
    bar.set_label(r"$h$ (nm)", fontsize=fs_labels)
    # Add title/ description
    ax.set_title(title, fontsize=fs_labels)
    # Axes ticks
    ax.set_xticks(np.linspace(0, np.max(thickness), 6))
    ax.tick_params(axis='x', labelsize=fs_ticks)
    # Saving figure
    if save:
        fig.savefig(f"../img/thickness_colour_{name}.svg", bbox_inches='tight')
    plt.show()

def plot_source(wavelengths, source_sd, name="", save=0):
    """ Plot scaled source spectral distribution """
    scale_sd = source_sd / np.max(source_sd)

    fs_labels = 12  # Fontsize of axes labels
    fs_ticks = 11  # Fontsize of tick labels
    lw = 1.4  # Linewidth
    width = 8.6 / 2.54  # Figure width (inches)
    fig, ax = plt.subplots(figsize=(1.5 * width, width), constrained_layout=True,
                           dpi=123.6)
    ax.plot(wavelengths, scale_sd, 'k-', lw=lw)
    ax.set_xlabel(r"$\lambda_{\mathrm{a}}$ (nm)", fontsize=fs_labels)
    ax.set_ylabel(r"Rel. $I_{\lambda \mathrm{s}}$", fontsize=fs_labels)
    # Set size of tick labels
    ax.tick_params(axis='both', labelsize=fs_ticks)
    # Remove top and right spines
    ax.spines.right.set_visible(False)
    ax.spines.top.set_visible(False)
    # Set limit on y axis
    ax.set_ylim(ymax=1.01)
    # Saving figure
    if save:
        fig.savefig(f"../img/source_{name}.svg", bbox_inches='tight')
    plt.show()
