""" IlluminatedSoapFilm class for interference calculations """

import numpy as np 
import os
# Import colour-science package 
import colour

from decorators import timeit


def detected_infinite(R, delta):
    """ Compute ratio of detected irradiance to the source irrdiance for an infinite
    number of interfering monochromatic waves.
    """
    F = 4 * R / (1 - R)**2
    tmp = F * np.sin(delta / 2)**2
    Id_Is = tmp / (1 + tmp)  # Id / Is
    return Id_Is

# def convert_Id_to_XYZ(Id):
#     """ Convert detected spectral irradiance distributions to XYZ tristimulus 
#     values. 
#
#     Inputs:
#     Id = each row contains the detected spectral irradiance distribution at a 
#     given film thickness
#     """
#     # Maximum spectral luminous efficacy (lm/W)
#     k = 683.002      
#
#     nh, nw = np.shape(Id)  # Rows for each thickness and cols for each wavelength
#     film_color_XYZ = np.zeros((nh, 3)) 
#     for count in range(nh):
#         detected_sd = colour.SpectralDistribution(Id[count], shape)
#         with colour.utilities.suppress_warnings(python_warnings=True):
#             XYZ = colour.sd_to_XYZ(detected_sd, cmfs, k=k)
#         film_color_XYZ[count] = XYZ

# def detected_monochromatic(R, T, delta, N):
#     """ Compute the detected irradiance to the source irradiance through the 
#     relation derived for monochromatic waves.
#     """
#     def loop(N):
#         """ Loop approach to allow for calculation with any N """
#         # Non-interference terms first 
#         sum1 = 0
#         for j in range(2, N+1):
#             sum1 += R**(2*j - 4)
#         # Interference terms 
#         sum2 = 0
#         for j in range(2, N+1):
#             sum2 += R**(j-2) * np.cos((j-1) * delta)
#         # Double sum 
#         sum3 = 0 
#         for j in range(2, N):
#             sum_ = 0
#             for m in range(j+1, N+1):
#                 sum_ += R**(m-2) * np.cos((m-j) * delta)
#             sum3 += R**(j-2) * sum_ 
#
#         return R * (1 + T**2 * sum1 - 2 * T * sum2 + 2 * T**2 * sum3)
#
#     if check:
#         detected_loop_N4 = loop(4)
#         detected_loop_N5 = loop(5)
#
#         # From manually expanding equation
#         detected_N4 =  R * (1 - 2 * T * (np.cos(delta) + R * np.cos(2 * delta)
#                         + R**2 * np.cos(3 * delta))
#                         + T**2 * (1 + 2 * R * np.cos(delta) + 2 * R**2 * np.cos(2 * delta)
#                                  + R**2 + 2 * R**3 * np.cos(delta) + R**4 ))
#
#
#         detected_N5 = R * (1 - 2 * T * (np.cos(delta) + R * np.cos(2 * delta)
#                         + R**2 * np.cos(3 * delta) + R**3 * np.cos(4 * delta))
#                         + T**2 * (1 + 2 * R * np.cos(delta) + 2 * R**2 * np.cos(2 * delta)
#                                  + 2 * R**3 * np.cos(3 * delta) + R**2 + 
#                                  + 2 * R**3 * np.cos(delta) + R**4 + 2 * R**4 * np.cos(2 * delta)
#                                  + 2 * R**5 * np.cos(delta) + R**6))
#
#         print(f"N = 4: {detected_loop_N4 = :.8f} and {detected_N4 = :.8f}")        
#         print(f"N = 5: {detected_loop_N5 = :.8f} and {detected_N5 = :.8f}")        
#         np.testing.assert_allclose(detected_loop_N4, detected_N4, rtol=0.0, atol=1e-12)
#         np.testing.assert_allclose(detected_loop_N5, detected_N5, rtol=0.0, atol=1e-12)
#         print("Loop verified")
#
#     return loop(N)

class ColourSoapFilm:
    def __init__(self, theta_air, nair, nfilm, shape, source_sd):
        """ Properties of the air-soap-air setup
        theta_air = angle of incidence of light ray from light source (rads)
        nair = absolute index of refraction of air 
        nfilm = absolute index of refraction of film 
        """
        self.theta_air = theta_air
        self.nair = nair
        self.nfilm = nfilm
        self.shape = shape
        self.wavelengths = shape.wavelengths
        self.source_sd = source_sd

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

    @timeit
    def interference_inf(self, thickness, polarisation="random"):
        """ Perform interference calculations by using the interference derived 
        for monochromatic waves at discrete wavelengths in the source with an 
        infinite number of interfering waves from a single incident lightwave
        from the source.

        Inputs:
        source = spectral distribution of the source, numpy array with spectral 
                 irradiance specified at wavelengths = 360, 361, 362, ..., 830 nm
        thickness = thickness values to compute the inteference for, numpy array (nm)

        """
        # Calculate reflectivity and transmissivity
        R_perp, T_perp, R_parr, T_parr, theta_film = self.fresnel_calc()

        # Perform interference calculations by using the monochromatic relation
        # at each wavelength at each film thickness 
        self.nh = len(thickness)  # Number of thickness values
        self.nw = len(self.wavelengths)  # Number of discrete wavelengths
        # Check source has same length as the number of wavelengths
        if len(self.source_sd) != self.nw:
            raise Exception(f"{len(self.source_sd) = } and {self.nw = }: they should be the same")

        self.Id = np.zeros((self.nh, self.nw))
        for ind_h, h in enumerate(thickness):
            for ind_w, wavelength in enumerate(self.wavelengths):
                delta = 4 * np.pi * h * self.nfilm * np.cos(theta_film) \
                        / (wavelength * self.nair)  # (rads)

                # If polarisation is perpendicular to plane of incidence
                if polarisation == "perp":
                    self.Id[ind_h, ind_w] = detected_infinite(R_perp, delta) * self.source_sd[ind_w]

                # If polarisation is parallel to plane of incidence
                elif polarisation == "parr":
                    self.Id[ind_h, ind_w] = detected_infinite(R_parr, delta) * self.source_sd[ind_w]

                # If light is randomly polarised
                elif polarisation == "random":
                    Id_perp = detected_infinite(R_perp, delta)
                    Id_parr = detected_infinite(R_parr, delta)
                    self.Id[ind_h, ind_w] = 0.5 * (Id_perp + Id_parr) * self.source_sd[ind_w]

                else:
                    raise ValueError("polarisation must be \"perp\", \"parr\" or \"random\"")

    @timeit
    def interference_inf_vectorised(self, thickness, polarisation="random"):
        """ Perform interference calculations by using the interference derived 
        for monochromatic waves at discrete wavelengths in the source with an 
        infinite number of interfering waves from a single incident lightwave
        from the source. Fully vectorised.

        Inputs:
        source = spectral distribution of the source, numpy array with spectral 
                 irradiance specified at wavelengths = 360, 361, 362, ..., 830 nm
        thickness = thickness values to compute the inteference for, numpy array (nm)

        """
        # Calculate reflectivity and transmissivity
        R_perp, T_perp, R_parr, T_parr, theta_film = self.fresnel_calc()

        # Perform interference calculations by using the monochromatic relation
        # at each wavelength at each film thickness 
        self.nh = len(thickness)  # Number of thickness values
        self.nw = len(self.wavelengths)  # Number of discrete wavelengths
        # Check source has same length as the number of wavelengths
        if len(self.source_sd) != self.nw:
            raise Exception(f"{len(self.source_sd) = } and {self.nw = }: they should be the same")

        h_2d = np.reshape(thickness, (self.nh, 1))
        h_repeat = np.repeat(h_2d, self.nw, axis=1)
        phase_shift = 4 * np.pi * self.nfilm * np.cos(theta_film) * h_repeat \
                      / (self.wavelengths * self.nair)


        # Coefficient of finesse used in interference calculation
        finesse_perp = 4 * R_perp / (1 - R_perp)**2
        finesse_parr = 4 * R_parr / (1 - R_parr)**2
        
        tmp_perp = finesse_perp * np.sin(phase_shift/2)**2
        tmp_parr = finesse_parr * np.sin(phase_shift/2)**2

        # If polarisation is perpendicular to plane of incidence
        if polarisation == "perp":
            self.Id = tmp_perp / (1 + tmp_perp) * self.source_sd

        # If polarisation is parallel to plane of incidence
        elif polarisation == "parr":
            self.Id = tmp_parr / (1 + tmp_parr) * self.source_sd

        # If light is randomly polarised
        elif polarisation == "random":
            Id_perp = tmp_perp / (1 + tmp_perp)
            Id_parr = tmp_parr / (1 + tmp_parr)
            self.Id = 0.5 * (Id_perp + Id_parr) * self.source_sd

        else:
            raise ValueError("polarisation must be \"perp\", \"parr\" or \"random\"")


    def convert_Id_to_XYZ(self):
        """ Convert detected spectral irradiance distribution at each wavelength
        to XYZ tristimulus values. 
        """
        # Load colour-matching functions
        cmfs = colour.MSDS_CMFS['CIE 1931 2 Degree Standard Observer']
        # Maximum spectral luminous efficacy (lm/W) - not necessary since scaling
        # by factor alpha in conversion to sRGB colourspace. It will be used so it
        # is explicit what k value is being used for testing purposes.
        k = 683.002      

        self.film_color_XYZ = np.zeros((self.nh, 3))  # cols are X, Y, Z
        for count in range(self.nh):
            detected_sd = colour.SpectralDistribution(self.Id[count], self.shape)
            with colour.utilities.suppress_warnings(python_warnings=True):
                XYZ = colour.sd_to_XYZ(detected_sd, cmfs, k=k)
            self.film_color_XYZ[count] = XYZ

        # return film_color_XYZ

    def convert_Id_to_XYZ_vectorised(self):
        """ Convert detected spectral irradiance distribution at each wavelength
        to XYZ tristimulus values using vectorised functions. 
        """
        # Loading in colour-matching functions
        file_dir = os.path.realpath(os.path.dirname(__file__))  # Dir of this file
        cmfs_file = os.path.join(file_dir, "cmfs/CIE_xyz_1931_2deg.csv")
        wave_cmfs = np.loadtxt(cmfs_file, delimiter=',')
        cmfs_values = wave_cmfs[:, 1:]

        # Maximum spectral luminous efficacy (lm/W) - not necessary since scaling
        # by factor alpha in conversion to sRGB colourspace. It will be used so it
        # is explicit what k value is being used for testing purposes.
        k = 683.002      

        # Calculating the *CIE XYZ* tristimulus values for each thickness 
        XYZ_p = self.Id[..., np.newaxis] * cmfs_values[np.newaxis, ...]  # 3rd dimension is X Y Z
        self.film_color_XYZ = k * np.sum(XYZ_p, axis=1)  # Sum along wavelengths

        # return film_color_XYZ

    def convert_XYZ_to_sRGB_vectorised(self, alpha):
        """ Convert XYZ tristimulus values to the sRGB colour space """
        # Conversion from XYZ to sRGB
        conversion_matrix_XYZ_to_sRGB = np.array([[3.2406, -1.5372, -0.4986], [-0.9689, 1.8758, 0.0415],
                                                [0.0557, -0.2040, 1.0570]])
        uncorrected_sRGB = np.transpose(np.matmul(conversion_matrix_XYZ_to_sRGB, 
                                        (self.film_color_XYZ / alpha).T))
        film_color_linear_sRGBclipped = np.clip(uncorrected_sRGB, 0.0, 1.0)
        # Apply gamma correction
        film_color_sRGBclipped = np.where(film_color_linear_sRGBclipped <= 0.0031308, 
                                          film_color_linear_sRGBclipped * 12.92,
                                          1.055 * film_color_linear_sRGBclipped ** (1 / 2.4) - 0.055)
        return film_color_sRGBclipped


        


# def interference_randomly_polarised(thickness, wavelengths, spectral_distribution, N):
#     """ Performs interference calculations for randomly polarised light. This 
#     implementation is non-vectorised.
#     """
#     Nw = len(wavelengths)  # Number of wavelengths considered 
#     for ind_h, h in enumerate(thickness):
#         for j in range(Nw):
#             delta = 4 * np.pi * h * nfilm * np.cos(theta_film) / (wavelengths[j] * nair)  # (rads)
#             detected_perp = detected_ratio(R_perp, T_perp, delta, N, check=0)  # Idperp/Isperp
#             detected_parr = detected_ratio(R_parr, T_parr, delta, N, check=0)  # Idparr/Isparr
#             detected_I = 0.5 * (detected_perp + detected_parr)  # Id/Is
#             detected_I_arr[ind_h, j] = detected_I 
