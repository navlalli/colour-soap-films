""" Calculate colours of a soap film illuminated by daylight """

import numpy as np
# Import colour-science package 
import colour 

# Import from this package
import interference
import plot_utils

def thickness_colour(save=0):
    """ Find the thickness-colour relationship for daylight: randomly polarised
    light with N = infinity.

    Mixture of vectorised and non-vectorised functions used.

    Inputs:
    save = boolean for deciding whether to save plot 

    """
    # Create soap film
    film = interference.ColourSoapFilm(theta_air, nair, nfilm, shape, source_sd)
    # Perform interference calculations
    film.interference_inf_vectorised(thickness=h, polarisation="random") 
    # Convert detected spectral irradiance distribution to XYZ colourspace for 
    # each thickness - would not be appropriate to use film.convert_Id_to_XYZ_vectorised()
    # here since the source is not specified at wavelengths = 360, 361, ..., 380 nm
    film.convert_Id_to_XYZ()
    # Scaling parameter so XYZ are in interval [0, 1]
    alpha = np.max(film.film_colour_XYZ) / 0.8 
    # Convert colours to sRGB colourspace
    film_colour_sRGBclipped = film.convert_XYZ_to_sRGB_vectorised(alpha)

    # Plot thickness-colour relationship
    plot_utils.plot_thickness_colour(h, film_colour_sRGBclipped, "daylight", save=save)

if __name__ == "__main__":

    # Global constants
    theta_air = 45.0 * np.pi / 180  # (rads)
    # Define source - using D65 illuminant to represent daylight
    D65 = colour.SDS_ILLUMINANTS['D65']
    shape = D65.shape
    wavelengths = shape.wavelengths
    source_sd = D65.values
    nw = len(source_sd)  # Number of discrete wavelengths in the source_sd
    # Absolute indices of refraction
    nair = np.full(nw, 1.003)
    nfilm = np.linspace(1.42, 1.40, nw)  # Arbitrary

    # Check out the source
    plot_utils.plot_source(shape.wavelengths, source_sd, "daylight", save=0)

    h = np.linspace(0.001, 1000, 1000)  # Thickness values of interest (nm)
 
    # Find variation of colour with thickness for daylight illumination
    thickness_colour(save=0)

