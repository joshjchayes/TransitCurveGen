'''
orbital_analysis.py

Functions to help calculate/estimate orbital parameters
'''
from platon.constants import M_sun, R_sun, R_jup
import numpy as np

def find_orbital_radius(T_planet, T_star, albedo):
    '''
    Calculates the orbital radius of the planet in units of stellar radii,
    assuming a circular orbit

    Parameters
    ----------
    T_planet : float
        The isothermal temperature of the planet's atmosphere. Note that
        we are assuming this is the same as the equilibrium temperatue.
    T_star : float
        The blackbody temeprature of the star.
    albedo : float
        The Bond albedo of the planet

    Returns
    -------
    orbital_radius : float
        The orbital radius of the planet in units of stellar radii,
        assuming a circular orbit.
    '''
    return (T_star ** 2 * np.sqrt(1 - albedo))/(2 * T_planet ** 2)


def find_orbital_period(orbital_radius, M_star, R_star, unit='seconds'):
    '''
    Calculates the orbital period of a circular orbit using Kepler's 3rd law

    Parameters
    ----------
    orbital_radius : float
        The orbital radius in units of stellar radii
    M_star : float
        The mass of the star in solar masses
    R_star : float
        The radius of the star in solar radii
    unit : str, optional
        The unit to return the period in. Accepted are
        - 'seconds'/'s'
        - 'hours'/'h'
        - 'days'/'d'
        - 'years'/'y'
    Default is seconds
    '''
    a = orbital_radius * R_star * R_sun

    G = 6.674e-11

    P = np.sqrt((4 * np.pi**2 * a**3)/(G * M_star * M_sun))

    unit = unit.lower()
    if unit in ['seconds', 's']:
        return P
    if unit in ['hours', 'h']:
        return P / 3600
    if unit in ['days', 'd']:
        return P / (3600 * 24)
    if unit in ['years', 'y']:
        return P / 31536000

    else:
        raise ValueError('Unrecognised unit {}'.format(unit))


def find_transit_duration(period, Rp, Rs, inc, orbital_radius):
    '''
    Calculates the transit duration

    Parameters
    ----------
    period : float
        The period of the planet orbit
    Rp : float
        The planet radius in jupiter radii
    Rs : float
        The star radus in solar radii
    orbital_radius :
        The orbital radius of the planet in units of host radii,
        assuming a circular orbit.

    Returns
    -------
    transit_duration : float
        The transit duration in the same unit as the period
    '''
    # Put things in SI units
    Rs = Rs * R_sun
    Rp = Rp * R_jup
    a = orbital_radius * Rs

    # Impact parameter
    b = a * np.cos(inc * (np.pi/180)) / Rs

    return (period/np.pi) * np.arcsin(np.sqrt((Rp + Rs)**2 - (b*Rs)**2)/a)
