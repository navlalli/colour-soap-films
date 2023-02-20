""" Calculate colours of a soap film illuminated by a blackbody of specified
temperature
"""

import numpy as np
# Import colour-science package 
import colour 

# Import from this package
import interference
import plot_utils

def thickness_colour(save=0):
    """ Find the thickness-colour relationship for an illuminated soap film 

    Inputs:
    save = boolean for deciding whether to save plot 

    """
    # Create soap film
    film = interference.ColourSoapFilm(theta_air, nair, nfilm, shape, source_sd)
    # Perform interference calculations
    film.interference_inf(thickness=h, polarisation="random") 
    # Convert detected spectral irradiance distribution to XYZ colourspace for 
    # each thickness
    film.convert_Id_to_XYZ()
    # Scaling parameter so XYZ are in interval [0, 1]
    alpha = np.max(film.film_colour_XYZ) / 0.8 
    # Conversion from XYZ to sRGB
    film_colour_sRGB = colour.XYZ_to_sRGB(film.film_colour_XYZ / alpha)
    # Clipped sRGB values outside [0, 1] to 0.0 or 1.0
    film_colour_sRGBclipped = np.clip(film_colour_sRGB, 0.0, 1.0)

    # Plot thickness-colour relationship
    plot_utils.plot_thickness_colour(h, film_colour_sRGBclipped, f"body{temp}K", save=save)

def thickness_colour_vectorised(save=0):
    """ Find the thickness-colour relationship for an illuminated soap film using
    only vectorised functions.

    Inputs:
    save = boolean for deciding whether to save plot 

    """
    # Create soap film
    film = interference.ColourSoapFilm(theta_air, nair, nfilm, shape, source_sd)
    # Perform interference calculations
    film.interference_inf_vectorised(thickness=h, polarisation="random") 
    # Convert detected spectral irradiance distribution to XYZ colourspace for 
    # each thickness
    film.convert_Id_to_XYZ_vectorised()
    # Scaling parameter so XYZ are in interval [0, 1]
    alpha = np.max(film.film_colour_XYZ) / 0.8 
    # Convert colours to sRGB colourspace
    film_colour_sRGBclipped = film.convert_XYZ_to_sRGB_vectorised(alpha)

    # Plot thickness-colour relationship
    plot_utils.plot_thickness_colour(h, film_colour_sRGBclipped, f"body{temp}K", save=save)


if __name__ == "__main__":
    # Global constants
    theta_air = 37.5 * np.pi / 180  # (rads)
    # Spectral distribution of an ideal blackbody of specified temperature
    shape = colour.SpectralShape(360, 830, 1)
    temp = 6500  # (K)
    source_sd = colour.sd_blackbody(temp, shape).values
    nw = len(source_sd)  # Number of discrete wavelengths in the source_sd
    # Absolute indices of refraction
    nair = np.full(nw, 1.0)
    nfilm = np.full(nw, 1.4)

    # Thickness values of interest (nm)
    h = np.linspace(0.001, 1000, 2000)  

    # Check out the source
    plot_utils.plot_source(shape.wavelengths, source_sd, f"body{temp}K", save=0)

    # Non-vectorised
    thickness_colour(save=0)
    # Vectorised
    thickness_colour_vectorised(save=0)
