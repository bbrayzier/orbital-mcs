#!/usr/bin/env python

"""
SolarSystem.py: TBD.
"""

from Common import *

AstroList = []

Sol = AstroObject("Sol", SolarMass, Colour='y')
AstroList.append(Sol)

#Earth = AstroObject("Earth", EarthMass, Position=np.array([aphel_m, 0., 0.]), Velocity=np.array([0., aphel_vel, 0.]), Colour='b') # Aphelion
Earth = AstroObject("Earth", EarthMass, Position=np.array([aphel_m, 0., 0.]), Colour='b') # Aphelion


Earth.orbit = EllipticalOrbit(a=1.00000011*AstroUnit, e=0.01671022, i=0.00005, omega=-11.26064, w=102.94719, v=100.46435)
Earth.orbit.UpdatePosition()
Earth.velocity = Earth.orbit.ProvideVelocity(SolarMass)
AstroList.append(Earth)


#Mars = AstroObject("Mars", MarsMass, Position=np.array([aphel_m, 0., 0.]), Velocity=np.array([0., aphel_vel, 0.]))
#AstroList.append(Mars)
