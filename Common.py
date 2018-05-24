#!/usr/bin/env python

"""
Common.py: TBD.
"""

__author__      = "Ben Brayzier"
__copyright__   = "Copyright 2018, Ben Brayzier"

### IMPORTS ###

import numpy as np
import random as rand
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from Constants import *

### MATH FUNCTIONS ###
def NormaliseVec(vector):
	"""Simple function to return a normalised vector"""
	vector.dtype = float
	norm = np.linalg.norm(vector)
	normalised = vector/norm
	return normalised

### PHYSICS FUNCTIONS ###
def GravForce(this_object, other_object):
	"""Takes two AstroObjects and returns the gravitational force vector between them"""
	separation = other_object.position - this_object.position
	mag_separation = np.linalg.norm(separation)
	mag_force = GravConst * this_object.mass * other_object.mass / mag_separation**2
	unit_vec = NormaliseVec(separation)
	force = mag_force * unit_vec
	return force

def GravAccel(this_object, other_object):
	"""Takes two AstroObjects and returns the gravitational acceleration vector between them"""
	separation = other_object.position - this_object.position
	mag_separation = np.linalg.norm(separation)
	mag_accel = GravConst * other_object.mass / mag_separation**2
	unit_vec = NormaliseVec(separation)
	accel = mag_accel * unit_vec
	return accel
	
### CLASSES ###
class AstroObject:
	"""
	AstroObject: A class to store an object in the simulation, whether that be star, planet, spacecraft or other body.
		Members: name, mass, position(= [0,0,0]), velocity(= [0,0,0]) and colour(= red)
	"""
	def __init__(self, Name, Mass, Position=np.zeros(3, dtype=float), Velocity=np.zeros(3, dtype=float), Colour='r'):
		self.name = Name
		self.mass = Mass
		self.position = Position
		self.velocity = Velocity
		self.colour = Colour

class EllipticalOrbit:
	"""
	EllipticalOrbit: A class to store the orbital parameters that describe an orbit, and methods that concern them.
		Members: semimajor axis, eccentricity, inclination, longitudinal ascending node, argument of periapsis and true anomaly
		Methods: UpdatePosition - 
		         ProvideVelocity
	"""
	def __init__ (self, a=1., e=0., i=0., omega=0., w=0., v=0.):
		self.semimajor = a
		self.eccentricity = e
		self.inclination = np.deg2rad(i)
		self.long_asc_node = np.deg2rad(omega)
		self.arg_periaps = np.deg2rad(w)
		self.true_anomaly = np.deg2rad(v)
		
	def RotatePlane(self, vector):
		"""
		Rotate a vector from the local elliptic plane to the reference plane described by Kepler elements
		"""
		
		cos_w = np.cos(self.arg_periaps)
		sin_w = np.sin(self.arg_periaps)
		cos_i = np.cos(self.inclination)
		sin_i = np.sin(self.inclination) 
		cos_omega = np.cos(self.long_asc_node)
		sin_omega = np.sin(self.long_asc_node) 

		arg_periaps_rotation = np.array([[cos_w, sin_w, 0.], [-sin_w, cos_w, 0.], [0., 0., 1.]])
		inclination_rotation = np.array([[1., 0., 0.], [0., cos_i, sin_i], [0., -sin_i, cos_i]])
		long_asc_node_rotation = np.array([[cos_omega, sin_omega, 0.], [-sin_omega, cos_omega, 0.], [0., 0., 1.]])

		full_rotation = np.dot(arg_periaps_rotation, np.dot(inclination_rotation, long_asc_node_rotation))

		vector = np.dot(full_rotation, vector) # Rotate ellipse relative to reference plane
		
		return vector
		
	def UpdatePosition(self, pos_type="saved"):
		# Determine eccentric anomaly
		if pos_type.lower() == "random":
			E = 2. * np.pi * rand.random()
		elif pos_type.lower() == "saved":
			E = 2. * np.arctan( 1 / np.sqrt((1 + self.eccentricity)/(1 - self.eccentricity)) * np.tan(self.true_anomaly / 2) )
		else:
			raise NameError("Incorrect position type argument: {}".format(pos_type))

		x = self.semimajor * np.cos(E) # x = a * cos(E)
		y = np.sqrt(self.semimajor**2 - self.eccentricity**2) * np.sin(E) # y = b * sin(E)
		position = np.array([x, y, 0.])  # (x, y, z) vector wrt own elliptic plane

		self.true_anomaly = np.arctan2(y, x)

#		cos_w = np.cos(self.arg_periaps)
#		sin_w = np.sin(self.arg_periaps)
#		cos_i = np.cos(self.inclination)
#		sin_i = np.sin(self.inclination) 
#		cos_omega = np.cos(self.long_asc_node)
#		sin_omega = np.sin(self.long_asc_node) 
#
#		arg_periaps_rotation = np.array([[cos_w, sin_w, 0.], [-sin_w, cos_w, 0.], [0., 0., 1.]])
#		inclination_rotation = np.array([[1., 0., 0.], [0., cos_i, sin_i], [0., -sin_i, cos_i]])
#		long_asc_node_rotation = np.array([[cos_omega, sin_omega, 0.], [-sin_omega, cos_omega, 0.], [0., 0., 1.]])
#
#		full_rotation = np.dot(arg_periaps_rotation, np.dot(inclination_rotation, long_asc_node_rotation))
#
#		position = np.dot(full_rotation, position) # Rotate ellipse relative to reference plane

		self.position = self.RotatePlane(position)


	def ProvideVelocity(self, mass_at_focus):
		# Insert vis-viva eqn. Hard bit is working out where the focus is i.e. r - the separation between body and focus
		focus_x = -self.semimajor * self.eccentricity # x = -c = -ae		
		focus_pos = np.array([focus_x, 0., 0.])
		print focus_pos
		focus_pos = self.RotatePlane(focus_pos)
		print focus_pos		
		# Calculate speed
		separation = np.linalg.norm(self.position - focus_pos)
		speed = np.sqrt(GravConst * mass_at_focus * (2 / separation - 1 / self.semimajor))
		print separation
		print speed
		
		
		# Eccentric anomaly to find velocity direction
		E = 2. * np.arctan( 1 / np.sqrt((1 + self.eccentricity)/(1 - self.eccentricity)) * np.tan(self.true_anomaly / 2) )
		vel_direction = np.array([-self.semimajor * np.sin(E), self.semimajor * np.sqrt(1 - self.eccentricity**2) * np.cos(E), 0.]) # [-a sin(E), b cos(E), 0]		
		
		velocity = self.RotatePlane(speed * (vel_direction / np.linalg.norm(vel_direction)))
		
		print velocity
		print np.array([0., aphel_vel, 0.])
		
		return velocity