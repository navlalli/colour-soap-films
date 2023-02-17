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

# Import speed of light
from scipy.constants import c 

def define_source(mean_wavelength, sigma_scale, plot=0):
    """ Create Gaussian source 

    Inputs:
    mean_wavelength = mean wavelength (nm)
    sigma_scale = sigma scale determines the source bandwidth

    """ 
    w = 2 * np.pi * c / nair / (wavelengths * 1e-9)  # Angular frequncies (rad/s)
    mean_w = 2 * np.pi * c / nair / (mean_wavelength * 1e-9)
    sigma = sigma_scale * mean_w
    source_sd = 1 / (sigma * (2 * np.pi)**0.5) \
                * np.exp(-0.5 * ((w - mean_w[0]) / sigma)**2)

    # Plot the spectral distribution of the source
    if plot:
        fig, ax = plt.subplots(constrained_layout=True)
        ax.plot(wavelengths, source_sd, 'k-', lw=1.4)
        ax.set_xlabel(r"$\lambda$ (nm)", fontsize=12)
        ax.set_ylabel(r"Spectral distribution", fontsize=12)
        plt.show()

    return source_sd

def gaussian_thickness_colour():
    """ Find the thickness-colour relationship for a soap film illuminated by 
    a Gaussian source with N = 5
    """
    # Create soap film
    film = interference.ColourSoapFilm(theta_air, nair, nfilm, shape, source_sd)
    # Perform interference calculations
    film.interference_N(thickness=h, N=5, polarisation="random") 
    # Convert detected spectral irradiance distribution to XYZ colourspace for 
    # each thickness
    film.convert_Id_to_XYZ_vectorised()
    # Scaling parameter so XYZ are in interval [0, 1]
    alpha = np.max(film.film_colour_XYZ) / 0.8 
    # Convert colours to sRGB colourspace
    film_colour_sRGBclipped = film.convert_XYZ_to_sRGB_vectorised(alpha)

    # Plot thickness-colour relationship
    plot_utils.plot_thickness_colour(h, film_colour_sRGBclipped, "Gaussian source")

if __name__ == "__main__":

    # Global constants
    theta_air = 37.5 * np.pi / 180  # (rads)
    nair = np.full(471, 1.003)
    nfilm = np.linspace(1.42, 1.40, 471)  # Arbitrary
    shape = colour.SpectralShape(360, 830, 1)
    wavelengths = shape.wavelengths
    source_sd = define_source(550, 0.2, plot=1)
    # Thickness values of interest (nm)
    h = np.linspace(0.01, 1500, 100)  

    # Find variation of colour with thickness with a Gaussian source
    gaussian_thickness_colour()


