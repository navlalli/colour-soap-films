""" Plotting utilities """ 

import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

def plot_thickness_colour(thickness, film_colour_sRGB, title=""):
    """ Plot showing the variation of colour with thickness """
    # Creates colourmap from RGB values
    custom_cmap = mpl.colors.ListedColormap(film_colour_sRGB)
    # sets min and max of colour bar labels
    norm = mpl.colors.Normalize(vmin=np.min(thickness), vmax=np.max(thickness))
    
    fs_label = 12
    fs_ticks = 11
    # fs_text = 11
    width = 12.9 / 2.54  # inches
    fig, ax = plt.subplots(figsize=(width, 0.30 * width), constrained_layout=True,
                            dpi=123.6)
    bar = mpl.colorbar.ColorbarBase(ax, cmap = custom_cmap, norm = norm, 
                                    orientation = 'horizontal')
    bar.set_label(r"$h$ (nm)", fontsize=fs_label)
    # Axes ticks size
    ax.tick_params(axis='x', labelsize=fs_ticks)
    # Title
    ax.set_title(title, fontsize=fs_label)

    # if save:
    #     fig.savefig(journal_dir + "quasicolourbar.eps", bbox_inches='tight')
    plt.show()
