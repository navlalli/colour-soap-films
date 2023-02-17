""" Test calculation of reflectivity and transmissivity """

import numpy as np
import colour
# Import from this package
import interference

def test_fresnel_calc():
    """ Test the implementation of Fresnels equations and calculation of reflectivity
    and transmissivity
    """ 
    theta_air = 30 * np.pi / 180
    nair = np.full(471, 1.0)
    nfilm = np.full(471, 1.4)
    source_shape = colour.SpectralShape(360, 830, 1.0)
    source_sd = np.ones(471)
    film = interference.ColourSoapFilm(theta_air, nair, nfilm, source_shape,
                                       source_sd)

    R_perp, T_perp, R_parr, T_parr, theta_film = film.fresnel_calc()

    # Reflection
    rdash_perp = (film.nfilm * np.cos(theta_film) - film.nair * np.cos(film.theta_air)) \
                 / (film.nfilm * np.cos(theta_film) + film.nair * np.cos(film.theta_air))
    Rdash_perp = rdash_perp**2

    # Transmission
    t_perp = 2 * film.nair * np.cos(film.theta_air) \
              / (film.nair * np.cos(film.theta_air) + film.nfilm * np.cos(theta_film))
    tdash_perp = 2 * film.nfilm * np.cos(theta_film) \
                / (film.nfilm * np.cos(theta_film) + film.nair * np.cos(film.theta_air))
    Tdash_perp  = (film.nair * np.cos(film.theta_air)) \
                  / (film.nfilm * np.cos(theta_film)) * tdash_perp**2

    # Check calculations for polarisation perpendicular to the plane of incidence
    # Check that R = Rdash
    np.testing.assert_allclose(R_perp, Rdash_perp, rtol=1e-5, atol=1e-8)
    # Check that T = Tdash = t * tdash
    np.testing.assert_allclose(T_perp, Tdash_perp, rtol=1e-5, atol=1e-8)
    ttdash_perp = t_perp * tdash_perp
    np.testing.assert_allclose(T_perp, ttdash_perp, rtol=1e-5, atol=1e-8)

    # Check that R + T = 1 and Rdash + Tdash = 1
    np.testing.assert_allclose(R_perp + T_perp, 1, rtol=1e-5, atol=1e-8)
    np.testing.assert_allclose(Rdash_perp + Tdash_perp, 1, rtol=1e-5, atol=1e-8)

    # Reflection
    rdash_parr = (film.nair * np.cos(theta_film) - film.nfilm * np.cos(film.theta_air)) \
                 / (film.nfilm * np.cos(film.theta_air) + film.nair * np.cos(theta_film))
    Rdash_parr = rdash_parr**2

    # Transmission
    t_parr = 2 * film.nair * np.cos(film.theta_air) \
             / (film.nair * np.cos(theta_film) + film.nfilm * np.cos(film.theta_air))
    tdash_parr = 2 * film.nfilm * np.cos(theta_film) \
                 / (film.nfilm * np.cos(film.theta_air) + film.nair * np.cos(theta_film))
    Tdash_parr  = (film.nair * np.cos(film.theta_air)) \
                  / (film.nfilm * np.cos(theta_film)) * tdash_parr**2

    # Check calculations for polarisation parallel to the plane of incidence
    # Check that R = Rdash
    np.testing.assert_allclose(R_parr, Rdash_parr, rtol=1e-5, atol=1e-8)
    # Check that T = Tdash = t * tdash
    np.testing.assert_allclose(T_parr, Tdash_parr, rtol=1e-5, atol=1e-8)
    ttdash_parr = t_parr * tdash_parr
    np.testing.assert_allclose(T_parr, ttdash_parr, rtol=1e-5, atol=1e-8)

    # Check that R + T = 1 and Rdash + Tdash = 1
    np.testing.assert_allclose(R_parr + T_parr, 1, rtol=1e-5, atol=1e-8)
    np.testing.assert_allclose(Rdash_parr + Tdash_parr, 1, rtol=1e-5, atol=1e-8)

