##########################################
## Modules for the 3D particle simulations
##########################################

import numpy as np
import os
from scipy.stats import norm
import time

import module_velocity_distr_3D ### functions to create different velocity distributions

#### create output folders
def func_dir_out(dir_add, nPart, nTime, dyN, speed):
	
	try:
		os.stat(dir_add)
	except:
		os.mkdir(dir_add)
	
	### subfolder :: positions of particles	- (periodic boundary conditions)
	dir_pos_data = (dir_add + 'positions_nPart_' + 
		str('{:0>8d}'.format(nPart)) +
		'_nTime' + str('{:0>4d}'.format(nTime)) + '/') 
	
	try:
		os.stat(dir_pos_data)
	except:
		os.mkdir(dir_pos_data)
	
	### subfolder :: continuous positions of particles 
	dir_posCon_data = (dir_add + 'positionsCon_nPart_' + 
		str('{:0>8d}'.format(nPart)) +
		'_nTime' + str('{:0>4d}'.format(nTime)) + '/') 
	try:
		os.stat(dir_posCon_data)
	except:
		os.mkdir(dir_posCon_data)

		
	dir_dynamics = (dir_add + 'dynamics_' + 
		str('{:0>8d}'.format(nPart)) +
		'_nTime' + str('{:0>4d}'.format(nTime)) + '/') 
	try:
		os.stat(dir_dynamics)
	except: 
		os.mkdir(dir_dynamics)
	

	### subfolder :: velocities of particles	
	dir_vel_data = (dir_add + 'velocites_Part_' + 
		str('{:0>8d}'.format(nPart)) +
		'_nTime' + str('{:0>4d}'.format(nTime)) + '/') 
	try:
		os.stat(dir_vel_data)
	except:
		os.mkdir(dir_vel_data)
	
	return dir_pos_data, dir_posCon_data, \
		dir_vel_data, dir_dynamics


##############################################################
## function to save all input parametes in text file
#############################################################
def save_sim_params(dir_add, BoxDim, nTime, dt, dyN, nPart, speed):
	file_sim_params = (dir_add + 'simulation_parameters.txt')
	str_file = ['BoxDimensions ' + \
		str(BoxDim[0]) + ' ' + str(BoxDim[1]) + ' ' + str(BoxDim[2]),\
		'nTime ' + str(nTime), \
		'dt ' + str(dt),\
		'dyN ' + str(dyN),\
		'nPart ' + str(nPart),\
		'speed ' + str(speed)]
	np.savetxt(file_sim_params, str_file , \
	fmt='%s', delimiter=' ', newline='\n')
	
	
	
###############################################################
## Select initial particle positions
##
## Input:
## BoxDim :: dimension of simulation box
##	BoxDim[0] = x-dim, BoxDim[1] = y-dim, BoxDim[2] = z-dim
## nPart :: number of particles

## Output: 
## xPos0, yPos0, zPos0 :: position coordinates
## phi0, theta0 :: setting direction for directed motion (spherical coor.)
################################################################
def initialize_particle_pos(BoxDim, nPart, dir_pos_data,dir_posCon_data,dir_vel_data):
	## random start position
	xPos0 = np.random.rand(nPart) * np.float(BoxDim[0]) ### initial x-Position
	yPos0 = np.random.rand(nPart) * np.float(BoxDim[1]) ### initial y-Position
	zPos0 = np.random.rand(nPart) * np.float(BoxDim[2]) ### initial z-Position

	# phi0 = 2. * np.pi * np.random.rand(nPart) ### azimuthal angle of initial particle displacement
	# theta0 = np.pi * np.random.rand(nPart) ### polar angle
	phi0 = 2. * np.pi * np.ones(nPart) ### azimuthal angle of initial particle displacement
	theta0 = np.pi * np.ones(nPart) ### polar angle

	###### save initial particle positions
	nt = 0 ## time step == 0
	savePos_Vel(xPos0,yPos0,zPos0,\
		0.,0.,0.,0., \
		nt,nPart, \
		dir_pos_data,dir_posCon_data,dir_vel_data)

	return xPos0, yPos0, zPos0, phi0, theta0

###############################################################
## Displace all particles 
##
## Input:
## BoxDim :: dimension of simulation box
##	BoxDim[0] = x-dim, BoxDim[1] = y-dim, BoxDim[2] = z-dim
## nPart :: number of particles
## speed :: mean speed of particles 
## dyN :: select type of implemented dynamics
## xPos0, yPos0, zPos0 :: position coordinates
## phi0, theta0 :: setting direction for directed motion (spherical coor.)
## dt :: step size of time steps (always = 1 ???)
## nt :: timestep counter
## v_distr :: gaussian distribution with sigma=1, for constant distribution of particle velocities
## dir_pos_data,dir_posCon_data,dir_vel_data :: directories for output
################################################################
def particle_displacement\
	(BoxDim, nPart, speed, dyN, \
	xPos, yPos, zPos, phi, theta, dt, nt, v_distr,\
	dir_pos_data,dir_posCon_data,dir_vel_data):


	##################### DYNAMICS ##########################################################	
	### different velocity distributions
	### constant velocity distribution	
	if dyN == 0:
		### directed motion 
		vx, vy, vz = \
		module_velocity_distr_3D.directed_motion(nPart,speed,phi,theta)
	elif dyN == 1:	
		### constant speed - not used anymore
		vx, vy, vz = \
		module_velocity_distr_3D.const_speed(nPart,speed)	
	elif dyN == 2:	
		### Maxwell-Boltzmann distribution -- leads to Brownian Motion
		vx, vy, vz = \
		module_velocity_distr_3D.brownian(nPart,speed)
	elif dyN == 3:
		### box distribution - not used anymore
		vx, vy, vz = module_velocity_distr_3D.box_distr(nPart,speed)
	elif dyN == 4:
		### just downwards motion in positive y-direction
		vx, vy, vz = \
		module_velocity_distr_3D.sedimentation(nPart, speed, phi, theta)
	elif dyN == 5:
		### downwards motion with velocity fluctuations 
		vx, vy, vz = \
		module_velocity_distr_3D.sedimentation_vDistr(nPart, speed, phi, theta)
	elif dyN == 6:
		### distribution of downwards motion with velocity fluctuations in x-dir 
		vx, vy, vz = \
		module_velocity_distr_3D.sedimentation_vyDistr_vxDiffusion(nPart, speed, phi, theta, v_distr)
		

	speed_test = np.sqrt(vx**2. + vy**2. + vz**2.) ## speed vector	
	m_speed_nt = np.mean(speed_test) ## mean spead at this time step
	
	##### update positions with periodic boundary conditions
	xPos = np.mod(xPos + dt * vx, BoxDim[0]) # periodic boundary conditions
	yPos = np.mod(yPos + dt * vy, BoxDim[1]) # periodic boundary
	
	#### z-direction - mirror position at plane z=0 or z=BoxDim[2]
	zPosMod = np.mod(zPos + dt * vz, BoxDim[2]) ## periodic boundary conditions
	zPosCon = zPos + dt * vz ## continuous motion crossing boundary
	##### zPosCross == 0 particle in box, zPosCross < 0 partilce crossing boundary
	zPosCross = zPosCon - zPosMod ## all particles crossing boundary have negative sign
	## = 1 if corssing, else 0
	boolCross = 1-np.float64(np.isclose(zPosCross, 0.0))   	
	### -1 crossing lower boundary, 0 still inside box, 1 crossing upper boundary	
	boolCrossUpDn = np.sign(zPosCross) * boolCross
	boolCrossUp = boolCross * (boolCrossUpDn + 1.) / 2. ## all crossing up bnd = 1
	boolCrossDn = boolCross * (boolCrossUpDn - 1.) / 2. ## all crossing lower = -1	
	##### update z-Pos: 	
	### case1 particle inside box :: zPos = zPosMod
	### case2 particle crossing upper boundary at z = BoxDim[2] :: zPos = zPosCross
	### case3 particle crossing lower boundary at z = 0 :: zPos = zPosCon * (-1)
	zPos = zPosMod * (1. - boolCross) + \
		zPosCross * boolCrossUp + \
		zPosCon * boolCrossDn 
	
	##### must be swopt 0 <--> pi for consecutive timesteps to maintain direction in dir mot.
	theta = np.abs(theta - np.pi * boolCross)


	# #### continuous positions for msd
	# xPosCon = xPosCon + dt * vx 
	# yPosCon = yPosCon + dt * vy


	######################################################################
	#### SAVE OUTPUT #####################################################
	savePos_Vel(xPos,yPos,zPos, \
		vx,vy,vz,speed, \
		nt,nPart, \
		dir_pos_data,dir_posCon_data,dir_vel_data)


	return xPos, yPos, zPos, phi, theta, m_speed_nt
	
#############################################################
## Implemented particle dynamics
## 
## Input: 
## BoxDim :: dimension of simulation box
##	BoxDim[0] = x-dim, BoxDim[1] = y-dim, BoxDim[2] = z-dim
## nPart :: number of particles
## nTime :: number of time steps
## dyN :: type of dynamics
## speed :: mean speed of particles
############################################################
def particle_dynamics\
	(BoxDim, nPart, nTime, dyN, speed, dt,\
	dir_pos_data, dir_posCon_data, dir_vel_data, dir_dynamics,\
	start_time):

	###################################################################################################
	#### generate initial particle positions
	xPos0,yPos0,zPos0,phi0,theta0 = \
		initialize_particle_pos(BoxDim, nPart, dir_pos_data,dir_posCon_data,dir_vel_data)
	
	xPos = xPos0
	yPos = yPos0
	zPos = zPos0
		
	phi = phi0
	theta = theta0

	# ##### continous positions for MSD calculations
	# xPosCon = xPos0
	# yPosCon = yPos0
	# zPosCon = zPos0
	
	
	###################################################################################################
	### allocate for dynamics output
	m_speed = np.zeros(nTime)

	### create gaussian distribution for velocity distribution of particles (ensemble NOT time)
	v_distr = norm.rvs(size=1*nPart, scale=1.) ## std = 1, mean = 0
	
	

	#####################################################################
	#### loop over timesteps to move the particles 
	# particle positions are updated every timestep
	for nt in (np.arange(nTime)+1): # can not be parallized in a simple way
		#### time processed
		if nt % 500 == 0:
			print ('working on timestep ' + str(nt) +
				' out of ' + str(nTime+1))
			elapsed_time = time.time() - start_time
			el_min, el_sec = divmod(elapsed_time, 60)
			el_hrs, el_min = divmod(el_min, 60) 
			print ('elapsed time: %d:%02d:%02d' % (el_hrs,el_min,el_sec))
	
		#### particle displacement
		xPos, yPos, zPos, phi, theta, m_speed_nt = particle_displacement\
			(BoxDim, nPart, speed, dyN, \
			xPos, yPos, zPos, phi, theta, dt, nt, v_distr, \
			dir_pos_data, dir_posCon_data, dir_vel_data)

		m_speed[nt-1] = m_speed_nt
			
	##########################################################	
	## end loop over time steps

	
	##############################################################################################
	### save mean speed to file
	file_dy_out = (dir_dynamics + 'dynamics_' +
		str('{:0>8d}'.format(nPart)) +
		'_Ntimestep' + str('{:0>4d}'.format(nTime)) + '.dat')
	
	np.savetxt(file_dy_out, np.transpose((np.arange(nTime),m_speed, np.ones(nTime)*speed)), \
	fmt='%e', delimiter=' ', newline='\n')



################################################
#### save positions and velocities
def savePos_Vel(xPos,yPos,zPos, \
	vx,vy,vz,speed, \
	nt,nPart, \
	dir_pos_data,dir_posCon_data,dir_vel_data):

	#### save positions to data file (accurate)
	file_pos_out = (dir_pos_data + 'particle_positions_nPart' +
		str('{:0>8d}'.format(nPart)) +
		'_timestep' + str('{:0>8}'.format(nt)) + '.dat')
	
	np.savetxt(file_pos_out, np.transpose((xPos,yPos,zPos)), \
	fmt='%.10f', delimiter=' ', newline='\n')




		
	#### save velocities to data file
	# v_abs = vx**2. + vy**2. + vz**2.
	v_mean = np.sqrt(vx**2. + vy**2. + vz**2.)
	
	file_vel_out = (dir_vel_data + 'particle_velocities_nPart' +
		str('{:0>8d}'.format(nPart)) +
		'_timestep' + str('{:0>8}'.format(nt)) + '.dat')
	
	np.savetxt(file_vel_out, np.transpose((vx,vy,vz,v_mean)), \
	fmt='%.10f', delimiter=' ', newline='\n')

	






