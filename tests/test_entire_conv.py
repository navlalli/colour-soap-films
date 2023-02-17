""" Test that the colours obtained using the vectorised and non-vectorised
methods are the same.
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
source_sd = colour.sd_blackbody(4300, source_shape).values
# Create soap film

def test_conv():
    """ Compare colours of an illuminated soap film computed using the vectorised
    and non-vectorised methods
    """
    film = interference.ColourSoapFilm(theta_air, nair, nfilm, source_shape, source_sd)
    film.nh = 100
    h = np.random.rand(film.nh) * 1500  # Thickness values of interest (nm)

    # Non-vectorised method
    film.interference_inf(thickness=h, polarisation="random") 
    film.convert_Id_to_XYZ()
    alpha = np.max(film.film_colour_XYZ) / 0.8 
    film_colour_sRGB = colour.XYZ_to_sRGB(film.film_colour_XYZ / alpha)
    film_colour_sRGBclipped = np.clip(film_colour_sRGB, 0.0, 1.0)

    # Vectorised method
    film.interference_inf_vectorised(thickness=h, polarisation="random") 
    film.convert_Id_to_XYZ_vectorised()
    alpha = np.max(film.film_colour_XYZ) / 0.8 
    film_colour_sRGBclipped_vec = film.convert_XYZ_to_sRGB_vectorised(alpha)

    # Check that vectorised and non-vectorised methods return the same result
    # Larger tolerance used since compounding the effect of multiple functions 
    np.testing.assert_allclose(film_colour_sRGBclipped_vec, film_colour_sRGBclipped, rtol=1e-5, atol=1e-4)
