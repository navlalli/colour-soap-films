""" Calculate colours of a soap film illuminated by a blackbody of specified
temperature
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
# Import colour-science package 
import colour 

# Import from this package
import sys
sys.path.append("/home/nav/navScripts/colour-soap-films/src/")
import interference
import plot_utils

# def plot_thickness_colour(thickness, film_colour_sRGB, title=""):
#     """ Plot showing the variation of colour with thickness """
#     # Creates colourmap from RGB values
#     custom_cmap = mpl.colors.ListedColormap(film_colour_sRGB)
#     # sets min and max of colour bar labels
#     norm = mpl.colors.Normalize(vmin=np.min(thickness), vmax=np.max(thickness))
#     
#     fs_label = 12
#     fs_ticks = 11
#     # fs_text = 11
#     width = 12.9 / 2.54  # inches
#     fig, ax = plt.subplots(figsize=(width, 0.30 * width), constrained_layout=True,
#                             dpi=123.6)
#     bar = mpl.colorbar.ColorbarBase(ax, cmap = custom_cmap, norm = norm, 
#                                     orientation = 'horizontal')
#     bar.set_label(r"$h$ (nm)", fontsize=fs_label)
#     # Axes ticks size
#     ax.tick_params(axis='x', labelsize=fs_ticks)
#     # Title
#     ax.set_title(title, fontsize=fs_label)
#
#     # if save:
#     #     fig.savefig(journal_dir + "quasicolourbar.eps", bbox_inches='tight')
#     plt.show()

def stat(arr, name=""):
    """ Stast of input array """
    print(f"{name}: {np.min(arr) = }")
    print(f"{name}: {np.mean(arr) = }")
    print(f"{name}: {np.max(arr) = }")

def thickness_colour():
    """ Find the thickness-colour relationship for an illuminated soap film """
    # Create soap film
    film = interference.ColourSoapFilm(theta_air, nair, nfilm, shape, source_sd)
    # Perform interference calculations
    film.interference_inf(thickness=h, polarisation="random") 
    stat(film.Id, "Id non-vec")
    # Convert detected spectral irradiance distribution to XYZ colourspace for 
    # each thickness
    film.convert_Id_to_XYZ()
    stat(film.film_colour_XYZ, "XYZ non-vec")

    # Scaling parameter so XYZ are in interval [0, 1]
    alpha = np.max(film.film_colour_XYZ) / 0.8 
    # Conversion from XYZ to sRGB
    film_colour_sRGB = colour.XYZ_to_sRGB(film.film_colour_XYZ / alpha)
    # Clipped sRGB values outside [0, 1] to 0.0 or 1.0
    film_colour_sRGBclipped = np.clip(film_colour_sRGB, 0.0, 1.0)
    stat(film_colour_sRGBclipped, "sRGB non-vec")

    # Plot thickness-colour relationship
    plot_utils.plot_thickness_colour(h, film_colour_sRGBclipped, f"body{temp}K", save=0)

def thickness_colour_vectorised():
    """ Find the thickness-colour relationship for an illuminated soap film using
    only vectorised functions.
    """
    # Create soap film
    film = interference.ColourSoapFilm(theta_air, nair, nfilm, shape, source_sd)
    # Perform interference calculations
    film.interference_inf_vectorised(thickness=h, polarisation="random") 
    stat(film.Id, "Id vec")
    # Convert detected spectral irradiance distribution to XYZ colourspace for 
    # each thickness
    film.convert_Id_to_XYZ_vectorised()
    # stat(film.film_colour_XYZ, "XYZ vec")
    # Scaling parameter so XYZ are in interval [0, 1]
    alpha = np.max(film.film_colour_XYZ) / 0.8 
    # Convert colours to sRGB colourspace
    film_colour_sRGBclipped = film.convert_XYZ_to_sRGB_vectorised(alpha)
    stat(film_colour_sRGBclipped, "sRGB vec")

    # Plot thickness-colour relationship
    plot_utils.plot_thickness_colour(h, film_colour_sRGBclipped, f"body{temp}K", save=0)


if __name__ == "__main__":
    # Global constants
    theta_air = 37.5 * np.pi / 180  # (rads)
    nair = np.full(471, 1.0)
    # nfilm = 1.41
    nfilm = np.full(471, 1.4)
    # Spectral distribution of an ideal blackbody of temperature 6500 K
    shape = colour.SpectralShape(360, 830, 1)
    temp = 6500  # (K)
    source_sd = colour.sd_blackbody(temp, shape).values

    # Thickness values of interest (nm)
    h = np.linspace(0.001, 1000, 500)  

    # Check out the source
    plot_utils.plot_source(shape.wavelengths, source_sd, f"body{temp}K", save=0)

    # Non-vectorised
    # thickness_colour()
    # Vectorised
    thickness_colour_vectorised()
