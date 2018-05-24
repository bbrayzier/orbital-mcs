#!/usr/bin/env python

"""
OrbitalModel.py: TBD.
"""

__author__      = "Ben Brayzier"
__copyright__   = "Copyright 2018, Ben Brayzier"


from Common import *
from SolarSystem import *


### MAIN ###

time_step = 0.02 * Day_s
no_steps = int(8 * Year_s / time_step)

def MidpointMethod(ObjList):
	no_bodies = len(ObjList)
	force_array = np.zeros((no_bodies, no_bodies, 3), dtype=float)
	pos_vec = np.zeros((no_bodies, 3), dtype=float)


	for i in range(no_bodies):
		for j in range(no_bodies):
			if i == j:
				continue
			elif i > j:
				force_array[i,j,:] = -1 * force_array[j,i,:]
			else:
				force_array[i,j,:] = GravForce(ObjList[i], ObjList[j])

		prev_velocity = ObjList[i].velocity
		prev_position = ObjList[i].position
		acceleration = np.sum(force_array[i], axis=0) / ObjList[i].mass # a = F / m
		ObjList[i].velocity = prev_velocity + acceleration * time_step # v(t+1) = v(t) + a(t+.5) * del_t
		ObjList[i].position = prev_position + (prev_velocity + ObjList[i].velocity) * time_step / 2. # r(t+1) = r(t) + ( v(t) + v(t+1) ) * del_t / 2
		pos_vec[i] = ObjList[i].position
		#print "Prev velocity: {}, New velocity: {}, New position: {} \n".format(prev_velocity, ObjList[i].velocity, ObjList[i].position)

	return ObjList, pos_vec

def update_plot(i, data, scat):
    scat.set_offsets(data[i])
    return scat

no_bodies = len(AstroList)
StoredPositions = np.zeros((no_steps, no_bodies, 3), dtype=float)

for t in range(no_steps):
	[ AstroList, StoredPositions[t] ] = MidpointMethod(AstroList)

#colours = np.random.random(no_bodies)
sizes = np.zeros(no_bodies, dtype=float)
colours = []
for i in range(no_bodies):
	sizes[i] = 100. * (AstroList[i].mass/SolarMass)**(1./3.)
	colours.append(AstroList[i].colour)

fig = plt.figure()
ax = plt.axes(xlim=(-1.25*AstroUnit, 1.25*AstroUnit), ylim=(-1.25*AstroUnit, 1.25*AstroUnit))
ax.set_aspect('equal', 'box')

scat = ax.scatter(StoredPositions[0,:,0], StoredPositions[0,:,1], c=colours, s=sizes)
x_axis = range(no_steps)
ani = animation.FuncAnimation(fig, update_plot, interval=25, frames=x_axis[::10], fargs=(StoredPositions[::10,:,:2], scat))
plt.show()
		
