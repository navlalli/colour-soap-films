""" ColourSoapFilm class for interference and colour calculations """

import numpy as np 
import os
# Import colour-science package 
import colour

# Import decorator for timing functions
from decorators import timeit

def detected_infinite(R, delta):
    """ Compute ratio of detected irradiance to the source irradiance for an infinite
    number of interfering monochromatic waves.

    Inputs:
    R = reflectivity of the film 
    delta = phase shift between each successive wave reaching the detector (rads)
    """
    F = 4 * R / (1 - R)**2
    tmp = F * np.sin(delta / 2)**2
    Id_Is = tmp / (1 + tmp)  # Id / Is

    return Id_Is

def detected_N(R, T, delta, N):
    """ Compute ratio of detected irradiance to the source irradiance for N
    interfering monochromatic waves.
    """
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

    Id_Is = R * (1 + T**2 * sum1 - 2 * T * sum2 + 2 * T**2 * sum3)

    return Id_Is

class ColourSoapFilm:
    def __init__(self, theta_air, nair, nfilm, shape, source_sd):
        """ Properties of the illuminated air-soap film-air setup
        theta_air = angle of incidence of a light ray from the light source (rads)
        nair = absolute index of refraction of air 
        nfilm = absolute index of refraction of the soap film 
        shape = colour.SpectralShape() from colour-science package
        source_sd = spectral distribution of the source
        """
        self.theta_air = theta_air
        self.nair = nair
        self.nfilm = nfilm
        self.shape = shape
        self.wavelengths = shape.wavelengths  # Wavelengths (nm)
        self.source_sd = source_sd
        self.nw = len(self.wavelengths)  # Number of discrete wavelengths
        # Check source has same length as the number of wavelengths
        if len(source_sd) != self.nw:
            raise Exception("source_sd should have same length as nw:"
                            f" {len(source_sd) = } and {self.nw = }")

    def fresnel_calc(self):
        """ Use Fresnel's formulas to return the reflectance and transmittance
        of the film 

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

        T_parr = (self.nfilm * np.cos(theta_film)) / (self.nair * np.cos(self.theta_air)) * t_parr**2

        return R_perp, T_perp, R_parr, T_parr, theta_film

    @timeit
    def interference_inf(self, thickness, polarisation="random"):
        """ Perform interference calculations using the interference equation
        derived for monochromatic waves at discrete wavelengths in the source.
        An infinite number of interfering waves from a single incident lightwave
        from the source is considered.

        Inputs:
        thickness = thickness values of interest, numpy array (nm)
        polarisation = polarisation of light from the source, str

        """
        # Calculate reflectivity and transmissivity
        R_perp, T_perp, R_parr, T_parr, theta_film = self.fresnel_calc()

        # Perform interference calculations by using the monochromatic relation
        # at each wavelength at each film thickness 
        self.nh = len(thickness)  # Number of thickness values
        self.Id = np.zeros((self.nh, self.nw))
        for ind_h, h in enumerate(thickness):
            # Calculate phase difference at a given film thickness
            delta = 4 * np.pi * h * self.nfilm * np.cos(theta_film) \
                    / (self.wavelengths * self.nair)  # (rads)

            for ind_w in range(self.nw):
                # If polarisation is perpendicular to plane of incidence
                if polarisation == "perp":
                    self.Id[ind_h, ind_w] = detected_infinite(
                                            R_perp[ind_w], delta[ind_w]) \
                                            * self.source_sd[ind_w]

                # If polarisation is parallel to plane of incidence
                elif polarisation == "parr":
                    self.Id[ind_h, ind_w] = detected_infinite(
                                            R_parr[ind_w], delta[ind_w]) \
                                            * self.source_sd[ind_w]

                # If light is randomly polarised
                elif polarisation == "random":
                    Id_perp = detected_infinite(R_perp[ind_w], delta[ind_w])
                    Id_parr = detected_infinite(R_parr[ind_w], delta[ind_w])
                    self.Id[ind_h, ind_w] = 0.5 * (Id_perp + Id_parr) \
                                            * self.source_sd[ind_w]

                else:
                    raise ValueError("polarisation must be \"perp\", \"parr\""
                                     " or \"random\"")

        return self.Id

    @timeit
    def interference_inf_vectorised(self, thickness, polarisation="random"):
        """ Perform interference calculations using the interference equation
        derived for monochromatic waves at discrete wavelengths in the source.
        An infinite number of interfering waves from a single incident lightwave
        from the source is considered and a fully vectorised approach is used.

        Inputs:
        thickness = thickness values of interest, numpy array (nm)
        polarisation = polarisation of light from the source, str

        """
        # Calculate reflectivity and transmissivity
        R_perp, T_perp, R_parr, T_parr, theta_film = self.fresnel_calc()

        # Perform interference calculations by using the monochromatic relation
        # at each wavelength at each film thickness 
        self.nh = len(thickness)  # Number of thickness values
        h_2d = np.reshape(thickness, (self.nh, 1))
        h_repeat = np.repeat(h_2d, self.nw, axis=1)
        delta = 4 * np.pi * h_repeat * self.nfilm * np.cos(theta_film) \
                / (self.wavelengths * self.nair)

        # Coefficient of finesse used in interference calculations
        finesse_perp = 4 * R_perp / (1 - R_perp)**2
        finesse_parr = 4 * R_parr / (1 - R_parr)**2
        
        tmp_perp = finesse_perp * np.sin(delta/2)**2
        tmp_parr = finesse_parr * np.sin(delta/2)**2

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

        return self.Id

    @timeit
    def interference_N(self, thickness, N, polarisation="random"):
        """ Perform interference calculations by using the interference equation
        derived for monochromatic waves at discrete wavelengths in the source 
        for N interfering waves from a single incident lightwave from the source.

        Inputs:
        thickness = thickness values of interest, numpy array (nm)
        N = number of interfering waves
        polarisation = polarisation of light from the source, str

        """
        # Calculate reflectivity and transmissivity
        R_perp, T_perp, R_parr, T_parr, theta_film = self.fresnel_calc()

        # Perform interference calculations by using the monochromatic relation
        # at each wavelength at each film thickness 
        self.nh = len(thickness)  # Number of thickness values
        self.Id = np.zeros((self.nh, self.nw))
        for ind_h, h in enumerate(thickness):
            # Calculate phase difference at a given film thickness
            delta = 4 * np.pi * h * self.nfilm * np.cos(theta_film) \
                    / (self.wavelengths * self.nair)  # (rads)

            for ind_w in range(self.nw):
                # If polarisation is perpendicular to plane of incidence
                if polarisation == "perp":
                    self.Id[ind_h, ind_w] = detected_N(
                                            R_perp[ind_w], T_perp[ind_w],
                                            delta[ind_w], N) * self.source_sd[ind_w]

                # If polarisation is parallel to plane of incidence
                elif polarisation == "parr":
                    self.Id[ind_h, ind_w] = detected_N(
                                            R_parr[ind_w], T_parr[ind_w],
                                            delta[ind_w], N) * self.source_sd[ind_w]

                # If light is randomly polarised
                elif polarisation == "random":
                    Id_perp = detected_N(R_perp[ind_w], T_perp[ind_w], delta[ind_w], N)
                    Id_parr = detected_N(R_parr[ind_w], T_parr[ind_w], delta[ind_w], N)
                    self.Id[ind_h, ind_w] = 0.5 * (Id_perp + Id_parr) * self.source_sd[ind_w]

                else:
                    raise ValueError("polarisation must be \"perp\", \"parr\" or \"random\"")

        return self.Id

    @timeit
    def convert_Id_to_XYZ(self):
        """ Convert detected spectral irradiance distribution at each thickness
        to XYZ tristimulus values. 
        """
        # Load colour-matching functions
        cmfs = colour.MSDS_CMFS['CIE 1931 2 Degree Standard Observer']
        # Maximum spectral luminous efficacy (lm/W) - not necessary since scaling
        # by factor alpha in conversion to sRGB colourspace. It is used for testing
        # purposes
        k = 683.002      

        self.film_colour_XYZ = np.zeros((self.nh, 3))  # cols are X, Y, Z
        for count in range(self.nh):
            detected_sd = colour.SpectralDistribution(self.Id[count], self.shape)
            with colour.utilities.suppress_warnings(python_warnings=True):
                XYZ = colour.sd_to_XYZ(detected_sd, cmfs, k=k)
            self.film_colour_XYZ[count] = XYZ

        return self.film_colour_XYZ

    @timeit
    def convert_Id_to_XYZ_vectorised(self):
        """ Convert detected spectral irradiance distribution at each thickness
        to XYZ tristimulus values using vectorised functions. 
        Requires the source_sd to be specified at wavelengths = 360, 361, ..., 830 nm
        as the cmfs used are at 1 nm intervals for 360 to 830 nm.
        """
        # Check source is specified at the correct wavelengths
        msg = "Specify source_sd at wavelengths = 360, 361, ..., 830 nm"
        np.testing.assert_equal(self.wavelengths, np.linspace(360, 830, 471),
                                err_msg=msg)

        # Loading in colour-matching functions
        file_dir = os.path.realpath(os.path.dirname(__file__))  # Dir of this file
        cmfs_file = os.path.join(file_dir, "cmfs/CIE_xyz_1931_2deg.csv")
        wave_cmfs = np.loadtxt(cmfs_file, delimiter=',')
        cmfs_values = wave_cmfs[:, 1:]
        # Maximum spectral luminous efficacy (lm/W) - not necessary since scaling
        # by factor alpha in conversion to sRGB colourspace. It is used for testing
        # purposes
        k = 683.002      

        # Calculating the *CIE XYZ* tristimulus values for each thickness 
        XYZ_p = self.Id[..., np.newaxis] * cmfs_values[np.newaxis, ...]
        self.film_colour_XYZ = k * np.sum(XYZ_p, axis=1)  # Sum along wavelengths

        return self.film_colour_XYZ

    @timeit
    def convert_XYZ_to_sRGB_vectorised(self, alpha=1.0):
        """ Convert XYZ tristimulus values to the sRGB colour space using vectorised
        functions

        Inputs:
        alpha = scaling factor to ensure X, Y and Z values are in the interval [0, 1]
        
        """
        # Conversion from XYZ to sRGB
        conversion_matrix_XYZ_to_sRGB = np.array([[3.2406, -1.5372, -0.4986],
                                                  [-0.9689, 1.8758, 0.0415],
                                                  [0.0557, -0.2040, 1.0570]])

        uncorrected_sRGB = np.transpose(np.matmul(conversion_matrix_XYZ_to_sRGB, 
                                        (self.film_colour_XYZ / alpha).T))

        # Clip values outside of interval [0, 1] to 0.0 or 1.0
        film_colour_linear_sRGBclipped = np.clip(uncorrected_sRGB, 0.0, 1.0)

        # Apply gamma correction
        film_colour_sRGBclipped = np.where(film_colour_linear_sRGBclipped <= 0.0031308, 
                                           film_colour_linear_sRGBclipped * 12.92,
                                           1.055 * film_colour_linear_sRGBclipped 
                                           ** (1 / 2.4) - 0.055)

        return film_colour_sRGBclipped
