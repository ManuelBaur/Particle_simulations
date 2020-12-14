##########################################
##  Module to generate different velocity distributions in 3D
##########################################


import numpy as np
from scipy.stats import norm

########################################################################################
#### constant speed with contant direction: sqrt(vx**2 + vy**2 + vz**2) = const. speed
#### Every particle fullfills directed motion with constant speed
#### The direction of motion is initialized for the first time step.

## INPUT:
## nPart :: number of particles
## speed :: of particles
## phi0 :: azimutahl angle of initial particle displacement
## theta0 :: polar angle 

def directed_motion(nPart,speed,phi,theta):	
	# rand_ang = 2. * np.pi * 0.25 * np.ones(nPart) # np.random.rand(nPart) # random angle in rad	
	### check definition of spherical coordinates
	vx = speed * np.sin(theta) * np.cos(phi)
	vy = speed * np.sin(theta) * np.sin(phi)
	vz = speed * np.cos(theta)

	return vx, vy, vz



########################################################################################
#### constant speed with random direction: sqrt(vx**2 + vy**2) = const. speed
## Every particle has a constant speed. Direction is changing every timestep.

## very artificial motion - not used anymore 

## INPUT:
## nPart :: number of particles
## speed :: 


def const_speed(nPart,speed):	
	phi = 2. * np.pi * np.random.rand(nPart) ### azimuthal angle of initial particle displacement
	theta = np.pi * np.random.rand(nPart) ### polar angle

	vx = speed * np.sin(theta) * np.cos(phi)
	vy = speed * np.sin(theta) * np.sin(phi)
	vz = speed * np.cos(theta)
	
	return vx, vy, vz


#########################################################################
### 	BROWNINAN MOTION
######################################################################## 
## help for norm.rvs:: 
###	https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.norm.html
## input must be nPart * 3 for all dimension.
## mean speed of 3D brownian motion:
## 	v_mean = sqrt(8*sigma^2 / pi)
def brownian(nPart,speed):
	## sigma selected to set speed == v_mean
	sigma = np.sqrt(speed**2. * np.pi / 8.)
	
	#### normal distribution	
	v3D = norm.rvs(size=3*nPart,scale=sigma)
	
	#### extract velocities in all three dimensions	
	vx = v3D[0:np.uint32(nPart)]
	vy = v3D[np.uint32(nPart):2*np.uint32(nPart)] # / 2.
	vz = v3D[2*np.uint32(nPart):3*np.uint32(nPart)]

	### mean speed
	v_abs = np.sqrt(vx**2. + vy**2. + vz**2.)
	v_mean = np.mean(v_abs)

	#### check particle velocity distribtution
	### plotting the values before and after the transformation
	# plt.figure()
	# plt.subplot(221) # the first row of graphs
	# plt.hist(vx,bins=np.int(np.sqrt(nPart)))     
	# plt.subplot(222)
	# plt.hist(vy,bins=np.int(np.sqrt(nPart)))
	# plt.subplot(223) # the second contains
	# plt.hist(vz,bins=np.int(np.sqrt(nPart)))     
	# plt.subplot(224)
	# plt.hist(v_abs,bins=np.int(np.sqrt(nPart)))
	# plt.show()
	# import pdb; pdb.set_trace()
	return vx, vy, vz	



#################################################
### Box distribution
### Funktion to displace particles by a random number (uniform distr)  

## INPUT:
## nPart :: number of particles
## partSpeedFac :: mean of the abs() of the velocities in the box distribution
def box_distr(nPart,partSpeedFac):
	## to reach right partSpeedFac
	factor = partSpeedFac * 2.0	
	
	### random number generator to initialize velocites
	vx = (np.random.rand(nPart) - 0.5) * factor 
	vy = (np.random.rand(nPart) - 0.5) * factor 
	vz = (np.random.rand(nPart) - 0.5) * factor 

	return vx, vy, vz



########################################################################################
#### constant speed in positive y-direction (downwards): 
#### There is NO velocity distribution.
#### sqrt(vx**2 + vy**2 + vz**2) = const. speed
## INPUT:
## nPart :: number of particles
## speed :: of particles
## phi0 :: azimutahl angle of initial particle displacement
## theta0 :: polar angle 

def sedimentation(nPart,speed,phi,theta):	
	#### set direction of displacement
	phi = 0.25 * 2. * np.pi * np.ones(nPart) ### azimuthal angle --> vx = 0 
	theta = 0.5 * np.pi * np.ones(nPart) ### polar angle --> vz = 0
	# rand_ang = 2. * np.pi * 0.25 * np.ones(nPart) # np.random.rand(nPart) # random angle in rad	
	### check definition of spherical coordinates
	vx = speed * np.sin(theta) * np.cos(phi)
	vy = speed * np.sin(theta) * np.sin(phi)
	vz = speed * np.cos(theta)

	return vx, vy, vz


########################################################################################
#### constant speed in positive y-direction (downwards): 
#### On top of that a diffuive term can be added in x- and y-direction
#### sqrt(vx**2 + vy**2 + vz**2) = const. speed
## INPUT:
## nPart :: number of particles
## speed :: of particles
## phi0 :: azimutahl angle of initial particle displacement
## theta0 :: polar angle 

def sedimentation_vDistr(nPart,speed,phi,theta):	
	#### set direction of displacement
	phi = 0.25 * 2. * np.pi * np.ones(nPart) ### azimuthal angle --> vx = 0 
	theta = 0.5 * np.pi * np.ones(nPart) ### polar angle --> vz = 0
	
	#### create brownian velocity distribution for fluctuations in x-direction
	[dvx, dvy, dvz] = brownian(nPart, speed)
	# import pdb; pdb.set_trace()

	# rand_ang = 2. * np.pi * 0.25 * np.ones(nPart) # np.random.rand(nPart) # random angle in rad	
	### check definition of spherical coordinates
	vx = speed * np.sin(theta) * np.cos(phi) + dvx
	vy = speed * np.sin(theta) * np.sin(phi) # + dvy / 2.
	vz = speed * np.cos(theta)

	return vx, vy, vz


########################################################################################
#### Velocity distribution on the particle ensemble in y-direction (downwards): 
#### This sedimenting velocity is constant over time for every particle.
#### On top of that a diffusion in x-direction can be switched on.
#### sqrt(vx**2 + vy**2 + vz**2) = const. speed
## INPUT:
## nPart :: number of particles
## speed :: of particles
## phi0 :: azimutahl angle of initial particle displacement
## theta0 :: polar angle 
## v_distr :: gaussian distribution with sig = 1 to generate time constant velocity distr
def sedimentation_vyDistr_vxDiffusion(nPart,speed,phi,theta,v_distr):	
	#### set direction of displacement
	phi = 0.25 * 2. * np.pi * np.ones(nPart) ### azimuthal angle --> vx = 0 
	theta = 0.5 * np.pi * np.ones(nPart) ### polar angle --> vz = 0
	
	#### create brownian velocity distribution for fluctuations in x-direction
	# dvx = norm.rvs(size=nPart, scale=1.25) ## updated every timestep
	[dvx, dvy, dvz] = brownian(nPart, speed*0.5)	
	dvy = dvy / 2.
	
	#### constant veloctiy distribution over time in y-direction
	dvy_cv20 = v_distr * speed * 0.20 ### every particle sediments with its own but constant velocity

	### spherical coordinates
	vx = speed * np.sin(theta) * np.cos(phi) + dvx ### diffusive in x-direction
	vy = speed * np.sin(theta) * np.sin(phi) + dvy + dvy_cv20 ### constant distribution
	vz = speed * np.cos(theta) + dvz

	

	return vx, vy, vz



