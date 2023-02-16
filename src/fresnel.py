""" Fresnels equations are used to calculate the reflectivity and transmissivity """

import numpy as np

def fresnel_calc(theta_air, nair, nfilm):
    """ Uses Fresnel's formulas for calculating the reflectance and transmittance

    Inputs:
    theta_air = angle of incidence of light ray from light source (rads)
    nair = absolute index of refraction of air 
    nfilm = absolute index of refraction of film 

    Outputs:
    R_perp = reflectivity for light polarised perpendicular to the plane of 
    incidence
    T_perp = transmissivity for light polarised perpendicular to the plane of 
    incidence
    R_parr = reflectivity for light polarised parallel to the plane of 
    incidence
    T_parr = transmissivity for light polarised parallel to the plane of 
    incidence
    theta_film = angle of ray inside film (rads)

    """
    # Snell's law 
    theta_film = np.arcsin(np.sin(theta_air) * nair / nfilm)  # (rads)
    theta_film_degs = theta_film * 180 / np.pi  # (degs)
    
    # Fresnel's formulas for perpendicular-polarised light 
    # Amplitude coefficients
    r_perp = (nair * np.cos(theta_air) - nfilm * np.cos(theta_film)) \
             / (nair * np.cos(theta_air) + nfilm * np.cos(theta_film))
    rdash_perp = (nfilm * np.cos(theta_film) - nair * np.cos(theta_air)) \
                 / (nfilm * np.cos(theta_film) + nair * np.cos(theta_air))
    
    t_perp = 2 * nair * np.cos(theta_air) \
              / (nair * np.cos(theta_air) + nfilm * np.cos(theta_film))
    tdash_perp = 2 * nfilm * np.cos(theta_film) \
                / (nfilm * np.cos(theta_film) + nair * np.cos(theta_air))

    # Reflectivity and transmissivity
    R_perp = r_perp**2
    Rdash_perp = rdash_perp**2

    T_perp = (nfilm * np.cos(theta_film)) / (nair * np.cos(theta_air)) * t_perp**2
    Tdash_perp  = (nair * np.cos(theta_air)) / (nfilm * np.cos(theta_film)) * tdash_perp**2

    # Checks 
    # Check 1: check that R = Rdash
    np.testing.assert_allclose(R_perp, Rdash_perp, rtol=0.0, atol=1e-8)

    # Check 2: check that T = Tdash = t tdash
    np.testing.assert_allclose(T_perp, Tdash_perp, rtol=0.0, atol=1e-8)
    ttdash_perp = t_perp * tdash_perp
    np.testing.assert_allclose(T_perp, ttdash_perp, rtol=0.0, atol=1e-8)

    # Check3: check that R + T = 1 and Rdash + Tdash = 1
    np.testing.assert_allclose(R_perp + T_perp, 1, rtol=0.0, atol=1e-8)
    np.testing.assert_allclose(Rdash_perp + Tdash_perp, 1, rtol=0.0, atol=1e-8)

    # Print results
    print("Perpendicular polarisation")
    print(f"\nr = {r_perp:.6f} r' = {rdash_perp:.6f} t = {t_perp:.6f} t' = {tdash_perp:.6f}")
    print(f"\nR = {R_perp:.6f} R' = {Rdash_perp:.6f} T = {T_perp:.6f} T' = {Tdash_perp:.6f}")

    # Fresnel's formulas for parallel-polarised light 
    # Amplitude coefficients
    r_parr = (nfilm * np.cos(theta_air) - nair * np.cos(theta_film)) \
              / (nair * np.cos(theta_film) + nfilm * np.cos(theta_air))
    rdash_parr = (nair * np.cos(theta_film) - nfilm * np.cos(theta_air)) \
                 / (nfilm * np.cos(theta_air) + nair * np.cos(theta_film))

    t_parr = 2 * nair * np.cos(theta_air) / (nair * np.cos(theta_film) + nfilm * np.cos(theta_air))
    tdash_parr = 2 * nfilm * np.cos(theta_film) / (nfilm * np.cos(theta_air) + nair * np.cos(theta_film))

    # Reflectivity and transmissivity
    R_parr = r_parr**2
    Rdash_parr = rdash_parr**2

    T_parr = (nfilm * np.cos(theta_film)) / (nair * np.cos(theta_air)) * t_parr**2
    Tdash_parr  = (nair * np.cos(theta_air)) / (nfilm * np.cos(theta_film)) * tdash_parr**2

    # Checks 
    # Check 1: check that R = Rdash
    np.testing.assert_allclose(R_parr, Rdash_parr, rtol=0.0, atol=1e-8)

    # Check 2: check that T = Tdash = t tdash
    np.testing.assert_allclose(T_parr, Tdash_parr, rtol=0.0, atol=1e-8)
    ttdash_parr = t_parr * tdash_parr
    np.testing.assert_allclose(T_parr, ttdash_parr, rtol=0.0, atol=1e-8)

    # Check3: check that R + T = 1 and Rdash + Tdash = 1
    np.testing.assert_allclose(R_parr + T_parr, 1, rtol=0.0, atol=1e-8)
    np.testing.assert_allclose(Rdash_parr + Tdash_parr, 1, rtol=0.0, atol=1e-8)

    # Print results
    print("Parallel polarisation")
    print(f"\nr = {r_parr:.6f} r' = {rdash_parr:.6f} t = {t_parr:.6f} t' = {tdash_parr:.6f}")
    print(f"\nR = {R_parr:.6f} R' = {Rdash_parr:.6f} T = {T_parr:.6f} T' = {Tdash_parr:.6f}")

    return R_perp, T_perp, R_parr, T_parr, theta_film

