'''
Spectrum objects for transitcurvegen

These objects are useful for keeping track of spectra, and can be combined
to give spectra with different spectral resolutions at different wavebands

'''
from .combined_spectrum import CombinedSpectrum
import matplotlib.pyplot as plt

class Spectrum:
    def __init__(self, wavelengths, depths, info_dict):
        '''
        The TransitCurveGen Spectrum objects are designed to hold information
        on a particular spectrum. They hold transit depth and wavelength info,
        along with the parameters used to generate them.
        Spectra with the same parameters can be combined to produce spectra
        spanning larger wavelength ranges with different spectral resolutions
        '''

        self.parameters = {}
        self.parameters['Rs'] = info_dict['Rs']
        self.parameters['Rp'] = info_dict['Rp']
        self.parameters['Mp'] = info_dict['Mp']
        self.parameters['T_planet'] = info_dict['T_planet']
        self.parameters['logZ'] = info_dict['logZ']
        self.parameters['CO_ratio'] = info_dict['CO_ratio']
        self.parameters['scattering_factor'] = info_dict['scattering_factor']
        self.parameters['cloudtop_pressure'] = info_dict['cloudtop_pressure']
        self.parameters['T_star'] = info_dict['T_star']
        self.parameters['Ms'] = info_dict['Ms']
        self.parameters['albedo'] = info_dict['albedo']

        self.wavelengths = wavelengths
        self.depths = depths

        self.min_wavelength = self.wavelengths.min()
        self.max_wavelength = self.wavelengths.max()

        self.resolution = info_dict['resolution']
        self.noise_level = info_dict['noise']

        self.spectra = self # This is here to help with CombinedSpectrum


    def combine(self, spectrum):
        '''
        Combines the Spectrum with another to produce a spectrum which covers
        multiple wavebands with different resolutions

        Parameters
        ----------
        spectrum : TransitCurveGen.Spectrum
            The spectrum with which to combine

        Returns
        -------
        combined_spectrum : TransitCurveGen.CombinedSpectrum
            The combined spectrum
        '''
        raise NotImplementedError()
        #return CombinedSpectrum(self, spectrum)

    def plot(self, savepath=None, show=True):
        '''
        Plots the spectrum nicely

        Parameters
        ----------
        savepath : str or None, optional
            If provided, will attempt to save the figure to this path. Default
            is None.
        show : bool, optional
            If True, will display the plot. Default is True.
        '''
        fig, ax = plt.subplots()

        if self.noise_level is None:
            # No error bars!
            ax.plot(self.wavelengths, self.depths, linewidth=2)

        else:
            err = np.ones(len(self.wavelengths)) * self.noise_level * 1e-6
            ax.errorbar(self.wavelengths, self.depths, err, fmt='x')

        ax.set_xscale('log')
        ax.set_xlabel('Wavelength / m')
        ax.set_ylabel('Transit depth $\Delta F$')

        if savepath is not None:
            try:
                fig.savefig(savepath)
            except Exception as e:
                print(e)
                print('Exception raised whilst saving figure. Figure is not saved')

        if show:
            plt.show()
