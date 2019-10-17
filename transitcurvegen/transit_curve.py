'''
transit_curve.py

Object to deal with Transit Curves

'''
import numpy as np
import matplotlib.pyplot as plt
import batman

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
        self.parameters['ecc'] = parameter_dict['ecc']

        # batman.TransitParams which can be used to generate the curve again
        self.batman_params = batman.TransitParams()
        self.batman_params.t0 = self.parameters['t0']
        self.batman_params.rp = self.parameters['Rp']
        self.batman_params.per = self.parameters['P']
        self.batman_params.a = self.parameters['a']
        self.batman_params.limb_dark = self.parameters['limb_darkening_model']
        self.batman_params.u = self.parameters['limb_darkening_params']
        self.batman_params.inc = self.parameters['inc']
        self.batman_params.w = self.parameters['w']
        self.batman_params.ecc = self.parameters['ecc']


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
        displacement = np.random.normal(0, noise_ppm*1e-6, size=len(self.flux))

        unc = np.ones(len(self.uncertainty)) * noise_ppm * 1e-6

        return TransitCurve(self.time, self.flux-displacement, unc, self.parameters)

    def plot(self, savepath=None, show=True, title=None, plot_true_curve=False):
        '''
        Displays the TransitCurve in a nice plot.

        Parameters
        ----------
        savepath : str or None, optional
            If provided, will attempt to save the figure to this path. Default
            is None.
        show : bool, optional
            If True, will display the plot. Default is True.
        title : str or None, optional
            If given, sets the plot title. Default is None
        plot_true_curve : bool, optional
            If True and there is noise added to the curve, will plot the actual
            noiseless light curve as well.
        '''

        fig, ax = plt.subplots()

        if self.uncertainty[0] == 0:
            # No uncertanties
            ax.plot(self.time, self.flux, label='Noiseless curve')

        else:
            ax.errorbar(self.time, self.flux, self.uncertainty, fmt='x',
                        label='Simulated data')

            if plot_true_curve:
                model = batman.TransitModel(self.batman_params, self.time)
                flux = model.light_curve(self.batman_params)

                ax.plot(self.time, flux, label='Noiseless curve')


        ax.legend()

        ax.set_xlabel('Time')
        ax.set_ylabel('Flux')
        if title is not None:
            ax.set_title(title)

        if savepath is not None:
            try:
                fig.savefig(savepath)
            except Exception as e:
                print(e)
                print('Exception raised whilst saving figure. Figure is not saved')

        if show:
            plt.show()
