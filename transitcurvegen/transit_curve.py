'''
transit_curve.py

Object to deal with Transit Curves

'''

class TransitCurve:
    def __init__(self, time, flux, uncertainty, parameter_dict):
        '''
        This is a TransitCurve - the flux, uncertainty and times, plus all the
        parameter information you might need. You can also use this object to
        add systematics.
        '''
        self.time = time
        self.flux = flux
        self.uncertainty = uncertainty

        self.parameters = {}

        self.parameters['Rp'] = parameter_dict['Rp']
        self.parameters['P'] = parameter_dict['P']
        self.parameters['t0'] = parameter_dict['t0']
        self.parameters['a'] = parameter_dict['a']
        self.parameters['inc'] = parameter_dict['inc']
        self.parameters['limb_darkening_model'] = parameter_dict['limb_darkening_model']
        self.parameters['limb_darkening_params']= parameter_dict['limb_darkening_params']
        self.parameters['w'] = parameter_dict['w']


    def add_gaussian_noise(self, noise_ppm):
        '''
        If uncertainties are 0, adds Gaussian noise at the ppm level specified

        Parameters
        ----------
        noise_ppm : float
            The noise in parts per million

        Returns
        -------
        transit_curve : TransitCurveGen.TransitCurve
            A new TransitCurve instance with noise added
        '''
        displacement = np.random.normal(0, noise_ppm*1e-6)

        unc = np.ones(len(self.uncertainty)) * noise_ppm * 1e-6

        return TransitCurve(self.times, self.flux-displacement, unc, self.parameters)
