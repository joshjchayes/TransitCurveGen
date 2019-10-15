'''
limb_darkening_handler

Module to deal with limb darkening and getting coefficients easily

'''

from ldtk import LDPSetCreator, BoxcarFilter
import numpy as np

class LimbDarkeningHandler:
    def __init__(self, spectrum, T_eff, logg, z):
        '''
        Object to generate realistic limb darkening parameters

        Parameters
        ----------
        spectrum : TransitCurveGen.Spectrum
            Spectrum we want the LD coeffs for
        T_eff : tuple
            The effective temperature of the host star, given as
            (value, uncertainty) pair
        logg : tuple
            The log_10 of the surface gravity of the host star, with gravity
            measured in cm/s2. Should be given as a (value, uncertainty) pair.
        z : tuple
            The metalicity of the host, given as a (value, uncertainty) pair.
        '''

        # Make the filters
        log_wl = np.log10(spectrum.wavelengths)
        delta  = log_wl[1] - log_wl[0]

        edges = [log_wl[0]-delta] + [wl + delta for wl in log_wl]
        edges = np.array(edges)

        self.filters = [BoxcarFilter('{}'.format(i), edges[i], edges[i+1]) for i in range(len(edges)-1)]

        print(len(self.filters), len(spectrum.wavelengths))


        sc = LDPSetCreator(teff=T_eff, logg=logg, z=z, filters=self.filters)

        self.ld_profiles = sc.create_profiles()

    def get_coeffs(self, method='quadratic', do_mc=True):
        '''
        Gets limb darkening coefficients for the given method

        Parameters
        ----------
        method : str, optional
            The limb darkening model to use. Allowed are linear, quadratic, and
            nonlinear. Default is quadratic.
        do_mc : bool, optional
            If true, will do some mcmc sampling to get good values

        Returns
        -------
        coeff : array_like, shape (n_filters, n_coeffs)
            The coefficients for each filter
        err : array_like, shape (n_filters, n_coeffs)
            The uncertainty on each of the coefficients
        '''
        if method == 'linear':
            coeff, err = self.ld_profiles.coeffs_ln(do_mc=do_mc)
        elif method == 'quadratic':
            coeff, err = self.ld_profiles.coeffs_qd(do_mc=do_mc)
        elif method == 'nonlinear':
            coeff, err = self.ld_profiles.coeffs_nl(do_mc=do_mc)
        else:
            raise ValueError('Unrecognised ld_model {}'.format(ld_model))

        return coeff, err
