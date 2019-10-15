'''
host_analysis

Module for functions needed to calculate parameters about the
'''
import numpy as np

def calculate_logg(host_mass, host_radius):
    '''
    Calculates log10(g) for a host in units usable by PyDTk. g is in cm s^-2

    Parameters
    ----------
    host_mass : float
        The host mass in solar masses
    host_radius : float
        The host radius in solar radii
    '''

    m_sun = 1.989e30
    r_sun = 6.957e8
    G = 6.674e-11

    return round(np.log10((host_mass * m_sun * G)/((host_radius * r_sun) ** 2) * 100), 2)
