""" Colours of a soap film when illuminated by the D65 illuminant """

import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
# Import colour-science package 
import colour 

# Import from this package
import sys
sys.path.append("/home/nav/navScripts/colour-soap-films/src/")
import interference

def plot_thickness_colour(thickness, film_color_sRGB, title=""):
    """ Plot showing the variation of colour with thickness """
    # Creates colormap from RGB values
    custom_cmap = mpl.colors.ListedColormap(film_color_sRGB)
    # sets min and max of color bar labels
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
    #     fig.savefig(journal_dir + "quasicolorbar.eps", bbox_inches='tight')
    plt.show()

def daylight_thickness_colour():
    """ Find the thickness-colour relationship for a soap film illuminated by 
    a daylight
    """
    # shape = colour.SpectralShape(360, 830, 1)
    # shape = colour.SpectralShape(360, 830, 1)
    # D65 = colour.SDS_ILLUMINANTS['D65']
    # print(f"{D65.shape = }")
    # print(f"{D65[550] = }") 
    # print(f"{D65.wavelengths[550] = }")
    # print(f"{D65.values[550] = }")
    # spectral_irradiance_source = light_source_sd.values  # W/sr/m^2/nm
    # wavelengths = D65.wavelengths
    # print(f"{D65.shape = }")
    # wavelengths = np.linspace(0, 1, 50)
    # spectral_values = D65.values

    # Approximating daylight by the spectral distribution of an ideal black body
    # of temperature 6500 K
    shape = colour.SpectralShape(360, 830, 1)
    temp = 6500  # (K)
    source_sd = colour.sd_blackbody(temp, shape)

    h = np.linspace(0.01, 1500, 10)  # Thickness values of interest (nm)

    theta_air = 37.5 * np.pi / 180  # (rads)
    nair = 1.0
    nfilm = 1.41
    # Create soap film
    film = interference.IlluminatedSoapFilm(theta_air, nair, nfilm, shape, source_sd)
    # Perform interference calculations
    film.interference_inf(thickness=h, polarisation="random") 
    # Convert detected spectral irradiance distribution to XYZ colourspace for 
    # each thickness
    film_color_XYZ = film.convert_Id_to_XYZ()

    # Scaling parameter so XYZ are in interval [0, 1]
    alpha = np.max(film_color_XYZ) / 0.8 

    # Conversion from XYZ to sRGB
    film_color_sRGB = colour.XYZ_to_sRGB(film_color_XYZ / alpha)
    # Clipped sRGB values outside [0, 1] to 0.0 or 1.0
    film_color_sRGBclipped = np.clip(film_color_sRGB, 0.0, 1.0)

    plot_thickness_colour(h, film_color_sRGBclipped, "Daylight illumination")


if __name__ == "__main__":

    daylight_thickness_colour()


