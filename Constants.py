import numpy as np

### CONSTANTS ###

# General # 
GravConst = 6.6740831e-11 # Units are m3kg-1s-2
AstroUnit = 149597870700 # Unitss are m

# Masses # 
EarthMass = 5.9722e24 # Units are kg
SolarMass = 1.98855e30 # Units are kg
#MarsMass = 

# Time #
Day_s = 24 * 3600 # One day in seconds
Year_s = 365.256363004 * Day_s # One year in seconds

# Earth's Orbit #
aphel_m = 1.0167 * AstroUnit
perihel_m = 0.98329 * AstroUnit
semimajor_m = 1.000001018 * AstroUnit

aphel_vel = np.sqrt( GravConst * SolarMass * ( 2 / aphel_m - 1 / semimajor_m ) ) # Vis-viva eqn.

# Mars' Orbit #
mars_aphel_m = 1.666 * AstroUnit
#mars_perihel_m

