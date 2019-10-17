'''
transit_curve_generator.py

An object to generate TransitCurves
'''
from .transit_curve import TransitCurve
from .orbital_analysis import find_orbital_period, find_orbital_radius, find_transit_duration
from .limb_darkening_handler import LimbDarkeningHandler
from .host_analysis import calculate_logg

import numpy as np
import batman

class TransitCurveGenerator:
    def __init__(self):
        '''
        The TransitCurveGenerator can be used to generate a set of TransitCurves
        from a transmission Spectrum, with capacity to add noise and systematics
        to each curve individually
        '''
        pass

    def curves_from_spectrum(self, spectrum, inc=90, ecc=0, w=90, T_eff=None,
                             logg=None, z=None, t0=None,
                             limb_darkening_model='quadratic',
                             limb_darkening_coeffs=None,
                             displace_time=True, noise_ppm=None):
        '''
        Generates a set of TransitCurves for the given spectrum

        Parameters
        ----------
        spectrum : TransitCurveGen.Spectrum
            The spectrum to generate the TransitCurves for.
        inc : float, optional
            The orbital inclination in degrees. Default is 90
        ecc : float, optional
            The eccentricity of the orbit. Note that period calculations are
            done assuming circular orbits, so change this at your own risk.
            Default is 0.
        w : float, optional
            The longitude of periastron in degrees. Default is 90
        T_eff : tuple or None, optional
            The effective temperature of the host star, given as
            (value, uncertainty) pair. If None, the T_star value from
            the spectrum will be used with a 1% error. Default is None.
        logg : tuple or None, optional
            The log_10 of the surface gravity of the host star, with gravity
            measured in cm/s2. Should be given as a (value, uncertainty) pair.
            If None, the host mass and radius will be used from the spectrum
            with a 5% error.
        z : tuple or None, optional
            The metalicity of the host, given as a (value, uncertainty) pair.
            If None, the metalicity of the atmosphere of the planet from the
            spectrum will be used with a 20% error.
        t0 : float or None, optional
            The time of t0 for the first transit. If provided, this is used to
            extrapolate t0 for the rest of the curves, simulating building up
            a spectrum through transit measurements over multiple epochs. If
            None, we assume that the transits were measured simultaneously
            through spectroscopic measurements.
        limb_darkening_model : str, optional
            The limb darkening model to use. Allowed are linear, quadratic, and
            nonlinear. Default is quadratic.
        limb_darkening_coeffs : array_like or None
            If provided, these will be used as limb darkening coefficients for
            all the transit curves. If None, PyLDTk will be used to generate
            realistic limb darkening parameter values for each waveband.
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
        light_curves : np.array, shape (n_depths,)
            The TransitCurves
        '''

        if spectrum.noise_level is not None:
            print('WARNING: Your spectrum has had simulated noise added to it. Be careful if you are wanting to use this to generate "true" light curves. We can add noise to the light curves, but you should do this only with a noiseless spectrum.')

        # Pull out some variables from the spectrum
        Ts = spectrum.parameters['T_star']
        Ms = spectrum.parameters['Ms']
        Rs = spectrum.parameters['Rs']
        Tp = spectrum.parameters['T_planet']


        # Deal with limb darkening
        if limb_darkening_coeffs is None:
            # Generate realistic limb darkening params
            if T_eff is None:
                T_eff = (Ts, Ts * 0.01)
            if logg is None:
                logg = calculate_logg(Ms, Rs)
                logg = (logg, logg * 0.05)
            if z is None:
                z = (10 ** spectrum.parameters['logZ'], 10 ** spectrum.parameters['logZ'] * 0.2)

            print(T_eff, logg, z)

            LDHandler = LimbDarkeningHandler(spectrum, T_eff, logg, z)

            ld_params, ld_errs = LDHandler.get_coeffs(method=limb_darkening_model)

        else:
            # Fixed LD params
            ld_params = [limb_darkening_coeffs for i in range(len(spectrum.wavelengths))]

        # Work out global variables (non-variant with wavelength or epoch)
        orbital_radius = find_orbital_radius(Tp, Ts, spectrum.parameters['albedo'])

        period = find_orbital_period(orbital_radius, Ms, Rs)

        # Sort out if the data is spectroscopic or not.
        if t0 is None:
            t0 = 0
            spectroscopic = True
        else:
            spectroscopic = False

        # Array to store the light curves
        light_curves = np.zeros(len(spectrum.depths), object)

        # Sort out the noise things
        if noise_ppm is not None:
            try:
                iterator = iter(noise)
            except:
                noise_ppm = np.ones(len(light_curves)) * noise_ppm

        # Loop through each data point in the spectrum and create a light curve

        for i, d in enumerate(spectrum.depths):
            print('Generating light curve {}/{}'.format(i+1, len(spectrum.wavelengths)))
            Rp = np.sqrt(d) * spectrum.parameters['Rs']

            # Make the batman parameters
            params = batman.TransitParams()
            params.t0 = t0
            params.per = period
            params.rp = Rp
            params.a = orbital_radius
            params.inc = inc
            params.ecc = ecc
            params.w = w
            params.limb_dark = limb_darkening_model
            params.u = ld_params[i]

            # What times are we evaluating the model at?
            # What is the expected transit duration?
            duration = find_transit_duration(period, Rp, Rs, inc, orbital_radius)

            # how far either side of t0 without displacement.
            delta_ti = duration/2 + 0.5 * duration
            delta_tf = duration/2 + 0.5 * duration

            if displace_time:
                # Now we move ∆ti and ∆tf randomly
                d1 = np.random.normal(0, delta_ti * 0.2)
                d2 = np.random.normal(0, delta_tf * 0.2)


                while (t0 + delta_ti + d1 >= t0 + delta_tf + d2) or t0-delta_ti+d1 >t0 or t0+delta_tf+d2 < t0:
                    #print(d1, d2)

                    d1 = np.random.normal(0, 0.01)
                    d2 = np.random.normal(0, 0.01)

                delta_ti += d1
                delta_tf += d2

            times = np.linspace(t0 - delta_ti, t0 + delta_tf, 200)

            # Calculate the transit light curve!
            model = batman.TransitModel(params, times)
            flux = model.light_curve(params)

            unc = np.zeros(len(times))

            parameter_dict = {'Rp':Rp, 'P':period, 't0':t0, 'a':orbital_radius,
                              'inc':inc, 'limb_darkening_params':ld_params[i],
                              'limb_darkening_model':limb_darkening_model,
                              'w':w, 'ecc':ecc}

            # Make the noiseless TransitCurve
            light_curves[i] = TransitCurve(times, flux, unc, parameter_dict)

            # Add Gaussian noise?
            if noise_ppm is not None:
                light_curves[i] = light_curves[i].add_gaussian_noise(noise_ppm[i])

            if not spectroscopic:
                t0 += period

        return light_curves
