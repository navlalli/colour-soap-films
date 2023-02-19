""" Calculate colours of a soap film illuminated by a Gaussian source """

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

def define_source(mean_wavelength, sigma_scale, name="", save=0):
    """ Create Gaussian source 

    Inputs:
    mean_wavelength = mean wavelength (nm)
    sigma_scale = sigma scale determines the source bandwidth
    name = string for saving source plot
    save = boolean, 0 for no saving and 1 for saving 

    Outputs:
    source_sd = spectral distribution of the source

    """ 
    # Define bandwidth relative to wavelength_fixed so it is the same for each 
    # Gaussian source used
    wavelength_fixed = 580  # (nm)
    sigma = sigma_scale * wavelength_fixed

    source_sd = 1 / (sigma * (2 * np.pi)**0.5) \
                * np.exp(-0.5 * ((wavelengths - mean_wavelength) / sigma)**2)

    print(f"{np.sum(source_sd) = }")
    # Check out the source
    plot_utils.plot_source(wavelengths, source_sd, f"{name}", save=save)

    return source_sd

def define_another_source(sigma_scale, name="", save=0):
    """ Create a source composed of the sum of Gaussians """
    # Define bandwidth relative to wavelength_fixed so it is the same for each 
    # Gaussian source used
    wavelength_fixed = 580  # (nm)
    sigma = sigma_scale * wavelength_fixed

    mean_wavelengths = np.array([440, 550])
    source_sd = np.zeros(len(wavelengths))
    for i in mean_wavelengths:
        source_sd += 1 / (sigma * (2 * np.pi)**0.5) \
                     * np.exp(-0.5 * ((wavelengths - i) / sigma)**2)

    # Check out the source
    plot_utils.plot_source(wavelengths, source_sd, f"{name}", save=save)

    return source_sd

def gaussian_thickness_colour(name="", save=0):
    """ Find the thickness-colour relationship for a soap film illuminated by 
    a specified light source emitting light polarised perpendicular to the plane
    of incidence. N = 5 used in the interference calculation.

    Mixture of vectorised and non-vectorised functions.

    Inputs;
    name = string for saving thickness-colour plot
    save = boolean, 0 for no saving and 1 for saving 
    """
    # Create soap film
    film = interference.ColourSoapFilm(theta_air, nair, nfilm, shape, source_sd)
    # Perform interference calculations
    film.interference_N(thickness=h, N=5, polarisation="perp") 
    # Convert detected spectral irradiance distribution to XYZ colourspace for 
    # each thickness
    film.convert_Id_to_XYZ_vectorised()
    # Scaling parameter so XYZ are in interval [0, 1]
    alpha = np.max(film.film_colour_XYZ) / 0.8 
    # Convert colours to sRGB colourspace
    film_colour_sRGBclipped = film.convert_XYZ_to_sRGB_vectorised(alpha)
    # Plot thickness-colour relationship
    plot_utils.plot_thickness_colour(h, film_colour_sRGBclipped, f"{name}", save=save)

if __name__ == "__main__":

    # Global constants
    theta_air = 35.0 * np.pi / 180  # (rads)
    shape = colour.SpectralShape(360, 830, 1)
    wavelengths = shape.wavelengths
    nw = len(wavelengths)
    nair = np.full(nw, 1.003)
    nfilm = np.linspace(1.42, 1.40, nw)  # Arbitrary
    mean_w = 660  # Mean wavelength (nm)
    sigma_s = 0.015  # Scaling factor determining frequency bandwidth
    source_sd = define_source(mean_w, sigma_s, f"mean{mean_w}nm{sigma_s}", save=0)

    h = np.linspace(0.001, 1000, 100)  # Thickness values of interest (nm)

    # Find variation of colour with thickness with a Gaussian source
    gaussian_thickness_colour(f"mean{mean_w}nm{sigma_s}", save=0)
    
    # Now considering source composed of the sum of Gaussians 
    source_sd = define_another_source(0.015, name="double_gaussian_440_550nm", save=0)
    # Find variation of colour with thicknes
    gaussian_thickness_colour("double_gaussian_440_550nm", save=0)
