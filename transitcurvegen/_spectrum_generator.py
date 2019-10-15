'''
SpectrumGenerator

Object to handle generation of spectra using PLATON
'''

import numpy as np
from platon.transit_depth_calculator import TransitDepthCalculator
from platon.constants import M_jup, R_jup, R_sun
from .spectrum import Spectrum

class SpectrumGenerator:
    def __init__(self):
        '''
        The SpectrumGenerator generates spectra at different resolutions and
        with noises
        '''
        self.calculator = TransitDepthCalculator()

    def generate_spectrum(self, Rs, Mp, Rp, T_planet, T_star, Ms, logZ=0,
                          CO_ratio=0.53, scattering_factor=1,
                          cloudtop_pressure=np.inf, resolution=None,
                          noise_level=None, max_wavelength=3e-5,
                          min_wavelength=3e-7, albedo=0.8):
        '''
        Generates a spectrum using the PLATON TransitDepthCalculator.

        Parameters
        ----------
        Rs : float
            Stellar radius in solar radii
        Mp : float
            Planet mass in Jupiter masses
        Rp : float
            Planet radius in Jupiter radii
        T_planet : float
            Temperature of the isothermal atmosphere in Kelvin
        T_star : float
            Blackbody temperature of the star
        Ms : float
            Mass of the star in solar masses
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
        noise_level : float or None, optional
            The noise on the data in ppm. If provided, Gaussian noise will be
            added to the spectrum. Default is None.
        max_wavelength : float, optional
            The maximum wavelength to consider in m. Default is 3e-5
        min_wavelength : float, optional
            The minimum wavelength to consider in m. Default is 3e-7

        Returns
        -------
        spectrum : TransitCurveGen.Spectrum
            A Spectrum object containing all the information on the spectrum
        '''

        # Generate the basic spectrum
        wavelengths, depths = self.calculator.compute_depths(Rs * R_sun,
                                                             Mp * M_jup,
                                                             Rp * R_jup,
                                                             T_planet,
                                                             logZ=logZ,
                                                             CO_ratio=CO_ratio,
                                                             scattering_factor=scattering_factor,
                                                             cloudtop_pressure=cloudtop_pressure)

        # Remove wavelengths not of interest
        # First cut off the large wavelengths from depths and wavelengths.
        # MUST be done in this order else depths can't reference wavelengths
        depths = depths[wavelengths <= max_wavelength]
        wavelengths = wavelengths[wavelengths <= max_wavelength]

        depths = depths[wavelengths >= min_wavelength]
        wavelengths = wavelengths[wavelengths >= min_wavelength]

        # Convolve to give a particular resolution
        if resolution is not None:
            wavelengths, depths = self._bin_spectrum(wavelengths, depths, resolution)

        # Add in noise if required
        if noise_level is not None:
            wavelengths, depths = self._add_noise(wavelengths, depths, noise_level)


        info_dict = {'Rs':Rs, 'Mp':Mp, 'Rp':Rp, 'T_planet':T_planet, 'logZ':logZ,
                     'CO_ratio':CO_ratio, 'scattering_factor':scattering_factor,
                     'cloudtop_pressure':cloudtop_pressure, 'resolution':resolution,
                     'noise':noise_level, 'T_star':T_star, 'Ms':Ms, 'albedo':albedo}

        return Spectrum(wavelengths, depths, info_dict)

    def _bin_spectrum(self, wavelengths, depths, R):
        '''
        Bins a given spectrum to a spectral resolution R

        Parameters
        ----------
        wavelengths : array_like, shape (n_wavelengths,)
            The wavelengths of each data point. We assume that these are
            logarithmically spaced.
        depths : array_like, shape (n_wavelengths,)
            The transit depth of each data point.
        R : int
            Spectral resolution to bin to.

        Returns
        -------
        wavelengths : np.array
            The centre wavelength of each new wavelength bin
        depths : np.array
            The transit depth in each of the new wavelength bins

        Notes
        -----
        The binning is done using convolution with a top-hat.
        '''
        # Take log base-10 of the wavelengths to work on resolution
        log_wl = np.log10(wavelengths)

        # Find the wavelength spacing in log space
        delta = round(log_wl[1] - log_wl[0], 7)

        # Find the original resolution of the spectrum
        R_original = 1 // delta

        if R > R_original:
            raise ValueError('Spectrum has resolution {}: cannot increase resolution to {}'.format(R_original, R))

        # Work out how many current wl bins we need for the given resolution
        nbins = int(1 // (R * delta))

        # Do the binning
        wavelengths = wavelengths[:(wavelengths.size // nbins) * nbins].reshape(-1, nbins).mean(axis=1)
        depths = depths[:(depths.size // nbins) * nbins].reshape(-1, nbins).mean(axis=1)

        return wavelengths, depths

    def _add_noise(self, wavelengths, depths, noise_ppm):
        '''
        Adds Gaussian noise to a given spectrum

        Parameters
        ----------
        wavelengths : array_like, shape (n_wavelengths,)
            The wavelengths of each data point. We assume that these are
            logarithmically spaced.
        depths : array_like, shape (n_wavelengths,)
            The transit depth of each data point.
        noise_ppm : int
            The noise in parts per million

        Returns
        -------
        wavelengths : np.array
            The wavelengths
        depths : np.array
            The depths with Gaussian noise added.
        '''

        disp_arr = np.zeros(len(depths))

        for i, d in enumerate(depths):
            disp_arr[i] = np.random.normal(0, noise_ppm * 1e-6)

        depths += disp_arr

        return wavelengths, depths
