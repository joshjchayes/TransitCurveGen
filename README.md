# TransitCurveGen

`TransitCurveGen` generates sets of realistic exoplanet transit curves at multiple wavelengths and epochs. When given atmospheric parameters, it first simulates a transmission spectrum before reverse engineering transit light curves and adding gaussian noise and (in future updates) basic systematics.

## Contents
[Requirements](#requirements)

[Installation](#installation)

[Basic usage](#basic_usage)

[Under the hood](#under_the_hood)

[Known issues](#known_issues)


<a name="requirements"></a>
## Requirements
Along with an installation of Python 3 (with the standard Conda distribution packages), TransitFit requires the following packages to be installed:

- [batman](https://www.cfa.harvard.edu/~lkreidberg/batman/index.html) - For light curve forward modelling

- [platon](https://github.com/ideasrule/platon) - For transmission spectrum modelling

- [Limb darkening toolkit (ldtk)](https://github.com/hpparvi/ldtk) - For realistic limb darkening parameters based on the host parameters


<a name="installation"></a>
## Installation
Currently, `TransitCurveGen` is not available through `pip`, though future versions will be. In order to install this package to your Python library, run
```
git clone https://github.com/joshjchayes/TransitCurveGen.git
cd TransitCurveGen
python setup.py install
```

<a name="basic_usage"></a>
## Basic usage
If you want to run the full TransitCurveGen generation, then you only need to use one function: ```generate_transit_curves```! This function requires quite a few arguments, which describe the atmosphere of the planet you are interested in simulating, some host parameters, and some information regarding the appearance of your data (noise etc).

The docstring for ```generate_transit_curves``` gives full information on how to use it, explaining the effect of each of the arguments, but to demonstrate just how easy it is to use, here is an example looking at a Jupiter-type planet with a 1000K atmosphere, orbiting a sun-like star. We will also assume that we have a spectral resolution of $R=30$ and can observe between 300 and 800 nm.

```python
import transitcurvegen as tcg

spectrum, light_curves = tcg.generate_transit_curves(1, 1, 5800, 1, 1, 1000, resolution=30, max_wavelength=8e-7, min_wavelength=3e-7)
 
```
This returns a ```Spectrum``` and a list of ```TransitCurves```, which are custom built objects, which you can find more about [here](#under_the_hood). In order to use this function, you must at least supply the radius, mass and temperature of the star and planet, in the order Rs, Ms, Ts, Rp, Mp, Tp. Planet mass/radius are in Jupiter masses/radii, whilst stellar mass/radius are in Solar masses/radii.

<a name="under_the_hood"></a>
## Under the hood
`TransitCurveGen` is built around two main types of object - ```Spectrum``` and ```TransitCurve```, along with a constructor object for each. 

The ```SpectrumGenerator``` takes planet and host parameters, resolution and wavelength limits and produces a `Spectrum` using the forward model provided by `PLATON`. There is also functionality to produce a `Spectrum` with added Gaussian noise, but this is not used in `TransitCurve` generation. 

The `TransitCurveGenerator` takes a `Spectrum` and calculates a `TransitCurve` for each data point in the spectrum. Since only mass, radii and temperatures of the host and planet are given, the `TransitCurveGenerator` must make some approximations to calculate some orbital parameters. The user can specify orbital inclination, but the radius of the orbit (assuming circular orbit), and the period are calculated to be consistent with the temperatures given, which is discussed [here](#realistic_orb_params). 

<a name="realistic_orb_params"></a>
### Calculating orbital parameters consistent with atmosphere and host temperatures

We assume that the isothermal atmospheric temperature is equal to the equilibrium temperature T_eq of the planet, and also assume that the planet's orbit is circular. Therefore 

<a href="https://www.codecogs.com/eqnedit.php?latex=T_{eq}&space;=&space;T_{*}\sqrt{\frac{R_*}{2a}}&space;\left(1-A\right)^{\frac{1}{4}}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?T_{eq}&space;=&space;T_{*}\sqrt{\frac{R_*}{2a}}&space;\left(1-A\right)^{\frac{1}{4}}" title="T_{eq} = T_{*}\sqrt{\frac{R_*}{2a}} \left(1-A\right)^{\frac{1}{4}}" /></a>

where a is the separation between host and planet, R_* is the host radius, A is the Bond albedo of the planet, and T_* is the blackbody temperature of the host.

This can be rearranged to give the orbital radius in units of the stellar radius:

<a href="https://www.codecogs.com/eqnedit.php?latex=\frac{a}{R_*}&space;=&space;\frac{T_*^2\sqrt{1-A}}{2T_{eq}^2}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\frac{a}{R_*}&space;=&space;\frac{T_*^2\sqrt{1-A}}{2T_{eq}^2}" title="\frac{a}{R_*} = \frac{T_*^2\sqrt{1-A}}{2T_{eq}^2}" /></a>

We can then use Kepler's 3rd law to calculate the period of this orbit 

<a href="https://www.codecogs.com/eqnedit.php?latex=P=&space;\sqrt{a^3\frac{4\pi^2}{GM_*}}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?P=&space;\sqrt{a^3\frac{4\pi^2}{GM_*}}" title="P= \sqrt{a^3\frac{4\pi^2}{GM_*}}" /></a>

which we then feed to `batman` in order to generate the `TransitCurves`.

<a name="requirements"></a>
## Known issues
Obviously, not everything always works. If you encounter an issue, first check the Issues page to see if it has been seen before. If not, raise an issue!

The main bug I am currently aware of is that realistic limb darkening with PyLDTk is not working. For now, I suggest that users stick to fixed limb darkening parameters.
