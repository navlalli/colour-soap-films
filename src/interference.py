""" IlluminatedSoapFilm class for interference calculations """

import numpy as np 

def detected_monochromatic(R, T, delta, N):
    """ Compute the detected irradiance to the source irradiance through the 
    relation derived for monochromatic waves.
    """
    def loop(N):
        """ Loop approach to allow for calculation with any N """
        # Non-interference terms first 
        sum1 = 0
        for j in range(2, N+1):
            sum1 += R**(2*j - 4)
        # Interference terms 
        sum2 = 0
        for j in range(2, N+1):
            sum2 += R**(j-2) * np.cos((j-1) * delta)
        # Double sum 
        sum3 = 0 
        for j in range(2, N):
            sum_ = 0
            for m in range(j+1, N+1):
                sum_ += R**(m-2) * np.cos((m-j) * delta)
            sum3 += R**(j-2) * sum_ 

        return R * (1 + T**2 * sum1 - 2 * T * sum2 + 2 * T**2 * sum3)

    if check:
        detected_loop_N4 = loop(4)
        detected_loop_N5 = loop(5)

        # From manually expanding equation
        detected_N4 =  R * (1 - 2 * T * (np.cos(delta) + R * np.cos(2 * delta)
                        + R**2 * np.cos(3 * delta))
                        + T**2 * (1 + 2 * R * np.cos(delta) + 2 * R**2 * np.cos(2 * delta)
                                 + R**2 + 2 * R**3 * np.cos(delta) + R**4 ))


        detected_N5 = R * (1 - 2 * T * (np.cos(delta) + R * np.cos(2 * delta)
                        + R**2 * np.cos(3 * delta) + R**3 * np.cos(4 * delta))
                        + T**2 * (1 + 2 * R * np.cos(delta) + 2 * R**2 * np.cos(2 * delta)
                                 + 2 * R**3 * np.cos(3 * delta) + R**2 + 
                                 + 2 * R**3 * np.cos(delta) + R**4 + 2 * R**4 * np.cos(2 * delta)
                                 + 2 * R**5 * np.cos(delta) + R**6))

        print(f"N = 4: {detected_loop_N4 = :.8f} and {detected_N4 = :.8f}")        
        print(f"N = 5: {detected_loop_N5 = :.8f} and {detected_N5 = :.8f}")        
        np.testing.assert_allclose(detected_loop_N4, detected_N4, rtol=0.0, atol=1e-12)
        np.testing.assert_allclose(detected_loop_N5, detected_N5, rtol=0.0, atol=1e-12)
        print("Loop verified")

    return loop(N)

class IlluminatedSoapFilm:
    def __init__(self, theta_air, nair, nfilm):
        """ Properties of the air-soap-air setup
        theta_air = angle of incidence of light ray from light source (rads)
        nair = absolute index of refraction of air 
        nfilm = absolute index of refraction of film 
        """
        self.theta_air = theta_air
        self.nair = nair
        self.nfilm = nfilm

    def fresnel_calc(self):
        """ Use Fresnel's formulas to return reflectance and transmittance

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
        theta_film = np.arcsin(np.sin(self.theta_air) * self.nair / self.nfilm)  # (rads)
        
        # Fresnel's formulas for perpendicular-polarised light 
        # Amplitude coefficients
        r_perp = (self.nair * np.cos(self.theta_air) - self.nfilm * np.cos(theta_film)) \
                 / (self.nair * np.cos(self.theta_air) + self.nfilm * np.cos(theta_film))
        
        t_perp = 2 * self.nair * np.cos(self.theta_air) \
                 / (self.nair * np.cos(self.theta_air) + self.nfilm * np.cos(theta_film))

        # Reflectivity and transmissivity
        R_perp = r_perp**2

        T_perp = (self.nfilm * np.cos(theta_film)) / (self.nair * np.cos(self.theta_air)) * t_perp**2

        # Fresnel's formulas for parallel-polarised light 
        # Amplitude coefficients
        r_parr = (self.nfilm * np.cos(self.theta_air) - self.nair * np.cos(theta_film)) \
                 / (self.nair * np.cos(theta_film) + self.nfilm * np.cos(self.theta_air))

        t_parr = 2 * self.nair * np.cos(self.theta_air) \
                 / (self.nair * np.cos(theta_film) + self.nfilm * np.cos(self.theta_air))

        # Reflectivity and transmissivity
        R_parr = r_parr**2

        T_parr = (self.nfilm * np.cos(theta_film)) \
                 / (self.nair * np.cos(self.theta_air)) * t_parr**2

        return R_perp, T_perp, R_parr, T_parr, theta_film




def interference_randomly_polarised(thickness, wavelengths, spectral_distribution, N):
    """ Performs interference calculations for randomly polarised light. This 
    implementation is non-vectorised.
    """
    Nw = len(wavelengths)  # Number of wavelengths considered 
    for ind_h, h in enumerate(thickness):
        for j in range(Nw):
            delta = 4 * np.pi * h * nfilm * np.cos(theta_film) / (wavelengths[j] * nair)  # (rads)
            detected_perp = detected_ratio(R_perp, T_perp, delta, N, check=0)  # Idperp/Isperp
            detected_parr = detected_ratio(R_parr, T_parr, delta, N, check=0)  # Idparr/Isparr
            detected_I = 0.5 * (detected_perp + detected_parr)  # Id/Is
            detected_I_arr[ind_h, j] = detected_I 
