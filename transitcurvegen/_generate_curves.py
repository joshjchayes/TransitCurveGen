'''
generate_curves

A wrapper function to generate transit curves from atmospheric info
'''

from ._transit_curve_generator import TransitCurveGenerator
from ._spectrum_generator import SpectrumGenerator
import numpy as np

def generate_transit_curves(Rs, Ms, Ts, Rp, Mp, Tp, logZ=0, CO_ratio=0.53,
                            scattering_factor=1, cloudtop_pressure=np.inf,
                            resolution=None, max_wavelength=3e-5,
                            min_wavelength=3e-7, albedo=0.8, inc=90,
                            ecc=0, w=90, limb_darkening_model='quadratic',
                            limb_darkening_coeffs=[0.1, 0.3],
                            displace_time=True, noise_ppm=None):
    '''
    Generates a set of transit curves. Creates a spectrum and works out transit
    curves from them.

    Parameters
    ----------
    Rs : float
        Stellar radius in solar radii
    Ms : float
        Mass of the star in solar masses
    Ts : float
        Blackbody temperature of the star
    Rp : float
        Planet radius in Jupiter radii
    Mp : float
        Planet mass in Jupiter masses
    Tp : float
        Temperature of the isothermal atmosphere in Kelvin
    logZ : float, optional
        base-10 logarithm of the metallicity in solar units. Default is 0.
    CO_ratio : float, optional
        C/O atomic ratio in the atmosphere. The default value is 0.53
        (solar) value.
    scattering_factor : float, optional
        Makes rayleigh scattering this many times as strong. Default is 1
    cloudtop_pressure : float, optional
        Pressure level (in PA) below which light cannot penetrate. Use
        np.inf for cloudless. Default is np.inf
    resolution : int or None, optional
        If provided, will produce spectra of the given spectral resolution.
        Default is None
    max_wavelength : float, optional
        The maximum wavelength to consider in m. Default is 3e-5
    min_wavelength : float, optional
        The minimum wavelength to consider in m. Default is 3e-7
    albedo : float,optional
        The Bond albedo of the planet. Default is 0.8
    inc : float, optional
        The orbital inclination in degrees. Default is 90
    ecc : float, optional
        The eccentricity of the orbit. Note that period calculations are
        done assuming circular orbits, so change this at your own risk.
        Default is 0.
    w : float, optional
        The longitude of periastron in degrees. Default is 90.
    limb_darkening_model : str, optional
        The limb darkening model to use. Allowed are linear, quadratic, and
        nonlinear. Default is quadratic.
    limb_darkening_coeffs : array_like or None
        If provided, these will be used as limb darkening coefficients for
        all the transit curves. If None, PyLDTk will be used to generate
        realistic limb darkening parameter values for each waveband. Default
        is [0.1, 0.3] for quadratic.
    displace_time : bool, optional
        If True, will shift calculations so that the center of the transit
        is not always at the centre of the TransitCurve times in order to
        parially simulate reality. Default is True.
    noise_ppm : int, array_like (shape (n_depths,)), or None, optional
        If not None, will add Gaussian noise to the light curves. If a
        single value is provided, this is the noise level in ppm that will
        be applied to each light curve. If an array is provided, it should
        be a list of noise levels in ppm which corresponds to each of the
        depths in the given spectrum. Default is None.

    Returns
    -------
    spectrum : TransitCurveGen.Spectrum
        The spectrum generated by the parameters given
    light_curves : np.array, shape (n_depths,)
        The TransitCurves which give each data point in the Spectrum
    '''

    spectrum_generator = SpectrumGenerator()
    transit_generator = TransitCurveGenerator()

    spectrum = spectrum_generator.generate_spectrum(Rs, Mp, Rp, Tp, Ts, Ms,
                                                    logZ, CO_ratio,
                                                    scattering_factor,
                                                    cloudtop_pressure,
                                                    resolution, None,
                                                    max_wavelength,
                                                    min_wavelength)

    light_curves = transit_generator.curves_from_spectrum(
        spectrum, inc, ecc, limb_darkening_model=limb_darkening_model,
        limb_darkening_coeffs=limb_darkening_coeffs,
        displace_time=displace_time, noise_ppm=noise_ppm)

    return spectrum, light_curves