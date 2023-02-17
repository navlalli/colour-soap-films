""" Test interference calculations """

import numpy as np
import colour
# Import from this package
import interference

theta_air = 35 * np.pi / 180
nair = np.full(471, 1.0)
nfilm = np.full(471, 1.4)
np.random.seed(50)
source_shape = colour.SpectralShape(360, 830, 1.0)
source_sd = np.random.rand(471) * 20e3
# Create soap film
film = interference.ColourSoapFilm(theta_air, nair, nfilm, source_shape, source_sd)

def test_detected_inf():
    """ Test interference calculation for N = infinity """ 
    h = np.random.rand(100) * 1500  # Thickness values of interest (nm)

    # Test for randomly polarised light
    Id = film.interference_inf(h, polarisation="random")
    Id_vec = film.interference_inf_vectorised(h, polarisation="random")
    # Check that vectorised and non-vectorised methods return the same result
    np.testing.assert_allclose(Id_vec, Id, rtol=1e-5, atol=1e-8)

    # Test for polarisation perpendicular to the plane of incidence
    Id = film.interference_inf(h, polarisation="perp")
    Id_vec = film.interference_inf_vectorised(h, polarisation="perp")
    # Check that vectorised and non-vectorised methods return the same result
    np.testing.assert_allclose(Id_vec, Id, rtol=1e-5, atol=1e-8)

    # Test for polarisation parallel to the plane of incidence
    Id = film.interference_inf(h, polarisation="parr")
    Id_vec = film.interference_inf_vectorised(h, polarisation="parr")
    # Check that vectorised and non-vectorised methods return the same result
    np.testing.assert_allclose(Id_vec, Id, rtol=1e-5, atol=1e-8)

def test_detected_N():
    """ Test interference calculations when specifying N """
    R_perp, T_perp, R_parr, T_parr, theta_film = film.fresnel_calc()

    delta = 0.875   # Arbitrary delta value

    N = 1
    # Perpendicular polarisation
    Id_Is_N1 = interference.detected_N(R_perp[50], T_perp[50], delta, N) 
    eq_N1 = R_perp  # Manually expanding sums in equation
    np.testing.assert_allclose(Id_Is_N1, eq_N1, rtol=1e-5, atol=1e-8)

    # Parallel polarisation
    Id_Is_N1 = interference.detected_N(R_parr[50], T_parr[50], delta, N) 
    eq_N1 = R_parr  # Manually expanding sums in equation
    np.testing.assert_allclose(Id_Is_N1, eq_N1, rtol=1e-5, atol=1e-8)

    N = 2
    def int_N2(R, T):
        """ Manually expanding sums in equation """
        return R * (1 - 2 * T * np.cos(delta) + T**2) 

    # Perpendicular polarisation
    Id_Is_N2 = interference.detected_N(R_perp[70], T_perp[70], delta, N) 
    eq_N2 = int_N2(R_perp[70], T_perp[70])
    np.testing.assert_allclose(Id_Is_N2, eq_N2, rtol=1e-5, atol=1e-8)

    # Parallel polarisation
    Id_Is_N2 = interference.detected_N(R_parr[70], T_parr[70], delta, N) 
    eq_N2 = int_N2(R_parr[70], T_parr[70])
    np.testing.assert_allclose(Id_Is_N2, eq_N2, rtol=1e-5, atol=1e-8)

    N = 3
    def int_N3(R, T):
        """ Manually expanding sums in equation """
        return R * (1 - 2 * T * (np.cos(delta) + R * np.cos(2 * delta))
                     + T**2 * (1 + 2 * R * np.cos(delta) + R**2)) 

    # Perpendicular polarisation
    Id_Is_N3 = interference.detected_N(R_perp[70], T_perp[70], delta, N) 
    eq_N3 = int_N3(R_perp[70], T_perp[70])
    np.testing.assert_allclose(Id_Is_N3, eq_N3, rtol=1e-5, atol=1e-8)

    # Parallel polarisation
    Id_Is_N3 = interference.detected_N(R_parr[70], T_parr[70], delta, N) 
    eq_N3 = int_N3(R_parr[70], T_parr[70])
    np.testing.assert_allclose(Id_Is_N3, eq_N3, rtol=1e-5, atol=1e-8)

    N = 4
    def int_N4(R, T):
        """ Manually expanding sums in equation """
        return R * (1 - 2 * T * (np.cos(delta) + R * np.cos(2 * delta)
                        + R**2 * np.cos(3 * delta))
                    + T**2 * (1 + 2 * R * np.cos(delta) + 2 * R**2 * np.cos(2 * delta)
                              + R**2 + 2 * R**3 * np.cos(delta) + R**4 ))

    # Perpendicular polarisation
    Id_Is_N4 = interference.detected_N(R_perp[250], T_perp[250], delta, N) 
    eq_N4 = int_N4(R_perp[250], T_perp[250])
    np.testing.assert_allclose(Id_Is_N4, eq_N4, rtol=1e-5, atol=1e-8)

    # Parallel polarisation
    Id_Is_N4 = interference.detected_N(R_parr[250], T_parr[250], delta, N) 
    eq_N4 = int_N4(R_parr[250], T_parr[250])
    np.testing.assert_allclose(Id_Is_N4, eq_N4, rtol=1e-5, atol=1e-8)

    N = 5
    def int_N5(R, T):
        """ Manually expanding sums in equation """
        return R * (1 - 2 * T * (np.cos(delta) + R * np.cos(2 * delta)
                                 + R**2 * np.cos(3 * delta) 
                                 + R**3 * np.cos(4 * delta))
                    + T**2 * (1 + 2 * R * np.cos(delta) 
                             + 2 * R**2 * np.cos(2 * delta)
                             + 2 * R**3 * np.cos(3 * delta) + R**2 + 
                             + 2 * R**3 * np.cos(delta) + R**4 
                             + 2 * R**4 * np.cos(2 * delta)
                             + 2 * R**5 * np.cos(delta) + R**6))

    # Perpendicular polarisation
    Id_Is_N5 = interference.detected_N(R_perp[150], T_perp[150], delta, N) 
    eq_N5 = int_N5(R_perp[150], T_perp[150])
    np.testing.assert_allclose(Id_Is_N5, eq_N5, rtol=1e-5, atol=1e-8)

    # Parallel polarisation
    Id_Is_N5 = interference.detected_N(R_parr[150], T_parr[150], delta, N) 
    eq_N5 = int_N5(R_parr[150], T_parr[150])
    np.testing.assert_allclose(Id_Is_N5, eq_N5, rtol=1e-5, atol=1e-8)
