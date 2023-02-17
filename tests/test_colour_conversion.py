""" Test conversion from detected spectral distribution to colour and conversion
from XYZ to sRGB colourspace
"""

import numpy as np
import colour
# Import from this package
import interference

theta_air = 35 * np.pi / 180
nair = np.full(471, 1.0)
nfilm = np.full(471, 1.4)
np.random.seed(22590)
source_shape = colour.SpectralShape(360, 830, 1.0)
source_sd = np.random.rand(471) * 20e3
# Create soap film
film = interference.ColourSoapFilm(theta_air, nair, nfilm, source_shape, source_sd)
film.nh = 100

def test_conversion_Id_to_XYZ():
    """ Test conversion from detected spectral distribution to colour """
    h = np.random.rand(film.nh) * 1500  # Thickness values of interest (nm)
    film.Id = np.random.rand(film.nh, film.nw) * 1000  # Create detected spectral distributions

    film_colour_XYZ = film.convert_Id_to_XYZ()
    film_colour_XYZ_vec = film.convert_Id_to_XYZ_vectorised()

    # Check that vectorised and non-vectorised methods return the same result
    np.testing.assert_allclose(film_colour_XYZ_vec, film_colour_XYZ, rtol=1e-5, atol=1e-8)

def test_conversion_XYZ_to_sRGB():
    """ Test conversion from XYZ to sRGB colourspace """
    film.film_colour_XYZ = np.random.rand(film.nh, 3)  # Create film colour in XYZ space

    # Conversion from XYZ to sRGB
    film_colour_sRGB = colour.XYZ_to_sRGB(film.film_colour_XYZ)
    film_colour_sRGB_clipped = np.clip(film_colour_sRGB, 0.0, 1.0)

    # Vectorised conversion from XYZ to sRGB
    film_colour_sRGB_clipped_vec = film.convert_XYZ_to_sRGB_vectorised(alpha=1.0)

    # Check that vectorised and non-vectorised methods return the same result
    np.testing.assert_allclose(film_colour_sRGB_clipped_vec, film_colour_sRGB_clipped, rtol=1e-5, atol=1e-8)



    



