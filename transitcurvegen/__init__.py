'''
transitcurvegen

A Python 3 package to enable creation of realistic transit light curves, based
upon transmission spectra simulated using PLATON

'''
__version__ = '0.1.1'
__name__ = 'transitcurvegen'


from ._spectrum_generator import SpectrumGenerator
from ._transit_curve_generator import TransitCurveGenerator
from ._generate_curves import generate_transit_curves
