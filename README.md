# The colour of soap films

This project shows how the colours exhibited by an illuminated soap film can be found by applying an interference relation derived for monochromatic waves at discrete wavelengths in the source. [Coherence effects on the colour of soap films]() provides a discussion of this method and its expected range of validity. All variables and equations implemented in `src/interference.py` are detailed in [coherence effects on the colour of soap films]().

Any of the example codes in `src` may be run with the following conda environment:
```
conda create -n soap-film numpy matplotlib
conda activate soap-film
pip install colour-science
```

## Examples

### Blackbody radiator
[`src/example_blackbody.py`](src/example_blackbody.py) shows how the colours of a soap film may be calculated when the light source is a blackbody radiator of specified temperature. This calculation considers an infinite number of interfering waves (see [paper]()). A fully vectorised implementation is provided in addition to a non-vectorised implementation. The non-vectorised method is more explicit, whilst the vectorised approach is **significantly** faster. The scaled spectral distribution of the source and the resulting thickness-colour relationships are shown for blackbodies of temperatures 3500 K, 6500 K and 9500 K:

| Source | Colour-thickness |
| :---: | :---: |
| 3500 K | 3500 K |
| ![Source for 3500 K](img/source_body3500K.svg) | ![Thickness-colour for 3500 K](img/thickness_colour_body3500K.svg) |
| 6500 K | 6500 K |
| ![Source for 6500 K](img/source_body6500K.svg) | ![Thickness-colour for 6500 K](img/thickness_colour_body6500K.svg) |
| 9500 K | 9500 K |
| ![Source for 9500 K](img/source_body9500K.svg) | ![Thickness-colour for 9500 K](img/thickness_colour_body9500K.svg) |

### Daylight 
[`src/example_daylight.py`](src/example_daylight.py) presents how the colours of a soap film illuminated by daylight may be calculated. This was done by using the D65 illuminant to represent daylight. This time, the coloured fringes are more widely spaced apart as a larger value was used for `theta_air`:

| Source | Colour-thickness |
| :---: | :---: |
| ![Source for daylight](img/source_daylight.svg) | ![Thickness-colour for daylight](./img/thickness_colour_daylight.svg) |

### Gaussian source 
[`src/example_gaussian_source.py`](src/example_gaussian_source.py) highlights how the colours of a soap film may be calculated when `source_sd` is specified by defining the numpy array through a calculation. In this case, a Gaussian source of chosen mean wavelength and frequency bandwidth emitting light polarised perpendicularly to the plane of incidence is used as the source. Furthermore, $N = 5$ was used in the interference calculations:

| Source | Colour-thickness |
| :---: | :---: |
| $\overline{\lambda} = 440$ nm | $\overline{\lambda} = 440$ nm |
| ![Gaussian source with mean wavelength at 440 nm](./img/source_mean440nm0.015.svg) | ![Thickness-colour for Gaussian source with mean wavelength at 440 nm](./img/thickness_colour_mean440nm0.015.svg) |
| $\overline{\lambda} = 550$ nm | $\overline{\lambda} = 550$ nm |
| ![Gaussian source with mean wavelength at 550 nm](./img/source_mean550nm0.015.svg) | ![Thickness-colour for Gaussian source with mean wavelength at 550 nm](./img/thickness_colour_mean550nm0.015.svg) |
| $\overline{\lambda} = 660$ nm | $\overline{\lambda} = 660$ nm |
| ![Gaussian source with mean wavelength at 660 nm](./img/source_mean660nm0.015.svg) | ![Thickness-colour for Gaussian source with mean wavelength at 660 nm](./img/thickness_colour_mean660nm0.015.svg) |

For a source with spectral distribution composed of the sum of two Gaussian distributions:

| Source | Colour-thickness |
| :---: | :---: |
| $\overline{\lambda} \in \lbrace 440, 550 \rbrace$ nm | $\overline{\lambda} \in \lbrace 440, 550 \rbrace$ nm |
| ![Sourced comprised of the sum of two Gaussians, centred at 440 and 550 nm](./img/source_double_gaussian_440_550nm.svg) | ![Sourced comprised of the sum of two Gaussians, centred at 440 and 550 nm](./img/thickness_colour_double_gaussian_440_550nm.svg) |


## License

This code is distributed under the [MIT License](LICENSE).
