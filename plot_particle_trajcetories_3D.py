################################################
##	particle trajectories in 3D	
#############################################

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.axes3d import Axes3D
import os

###################################################################################################
## parameters of fake images:
nPart = 1000 ## number of particles
dyN = 0 ## type of dynamics (1 = const speed, rand dirction, 2 = Max Boltz distr
partSpeedFac = 2 ## multiplication factor on mean speed of particles
noise = 1 ## Noise added to artificial images (1 = no noise, 2 = random px noise as in experiments)
nTime = 2000 ## number of time steps, size of steps = 1

### plot parameters
nTplot = 500 ## number of timesteps in plot
nPartPlot = 5 ## number of particles in plot

############################### path position data of particles 
f_meas = ('../../3D_ParameterScan1/' +
	'runs_nPart_' + str('{:d}'.format(nPart)) +
	'_dynamics_' + str('{:d}'.format(dyN)) +
	'_speed_' + str('{:d}'.format(partSpeedFac)) + '/')

dir_in = ('positions_nPart_' + 
	str('{:0>8d}'.format(nPart)) +
	'_nTime' + str('{:0>4d}'.format(nTime)) + '/') 
	

### directory for msd output - write to dynamics folder
dir_out = ('dynamics_' + str('{:0>8d}'.format(nPart)) + '_nTime' + str('{:0>4d}'.format(nTime)) + '/')

#### naming of input data files: exact with float number and rounded (floor) to pixel grid
nam_data_in = ('particle_positions_nPart' + str('{:0>8d}'.format(nPart)) + '_timestep') # exact data

try:
	os.stat(f_meas + dir_out)
except:
	os.mkdir(f_meas + dir_out)

####### allocate vector with all positions over time
##### down: positions over time, right: x and y pos, depth: number of particles
pos_vs_t = np.zeros((nTplot,3,nPart))

###############################
##### loop over all time steps
for nt in np.arange(nTplot): # nTime):
	
	
	### read data (exact and on px-grid)
	nam_data_t = (f_meas + dir_in + nam_data_in + str('{:0>4d}'.format(nt)) + '.dat')
	
	pos_t = np.loadtxt((nam_data_t))
	pos_vs_t[nt,0,:] = pos_t[:,0]
	pos_vs_t[nt,1,:] = pos_t[:,1]
	pos_vs_t[nt,2,:] = pos_t[:,2]
	
	




##################################################
#### settings for plot
################# plot and save data of msd
########## produce LaTex like fonts
font = {'family' : 'serif',
	'weight' : 'normal',
	'size' : 18}
plt.rc('text', usetex=True)
plt.rc('font', **font) # family='serif')
plt.rc
fig = plt.figure(figsize=(12,12))
ax = fig.add_subplot(221, projection='3d')

### loop over particles
for nnPart in range(nPartPlot): # range(nPart):
	ax.scatter(pos_vs_t[:,0,nnPart],pos_vs_t[:,1,nnPart],pos_vs_t[:,2,nnPart],\
		marker='.', s=1) 
	## green start point of trajectory
	ax.scatter(pos_vs_t[0,0,nnPart],pos_vs_t[0,1,nnPart],pos_vs_t[0,2,nnPart],\
		marker='o', s=10, color='C2') 	
	## red end point of trajectory
	ax.scatter(pos_vs_t[nTplot-1,0,nnPart],pos_vs_t[nTplot-1,1,nnPart],pos_vs_t[nTplot-1,2,nnPart],\
		marker='o', s=10, color='C3') 
	
ax.set_xlabel('$x$-position (Px)')
ax.set_ylabel('$y$-position (Px)')
ax.set_zlabel('$z$-position (Px)')
ax.axis('equal')
ax.set_aspect('equal','box')
ax.set_xticks(np.linspace(0,500,3))
ax.set_yticks(np.linspace(0,500,3))
ax.set_zticks(np.linspace(0,1000,5))
ax.set_xlim([0,512])
ax.set_ylim([0,512])
ax.set_zlim([0,1024])
# plt.show()



###############################################################################
######### projection to x-z plane
ax2 = fig.add_subplot(222)


### loop over particles
for nnPart in range(nPartPlot): # range(nPart):
	ax2.scatter(pos_vs_t[:,0,nnPart],pos_vs_t[:,2,nnPart],\
		marker='.', s=1) 
	## green start point of trajectory
	ax2.scatter(pos_vs_t[0,0,nnPart],pos_vs_t[0,2,nnPart],\
		marker='o', s=10, color='C2') 	
	## red end point of trajectory
	ax2.scatter(pos_vs_t[nTplot-1,0,nnPart],pos_vs_t[nTplot-1,2,nnPart],\
		marker='o', s=10, color='C3') 
	
ax2.set_xlabel('$x$-position (Px)')
ax2.set_ylabel('$z$-position (Px)')
ax2.axis('equal')
ax2.set_aspect('equal','box')
ax2.set_xticks(np.linspace(0,500,3))
ax2.set_yticks(np.linspace(0,1000,5))
ax2.set_xlim([0,512])
ax2.set_ylim([0,1024])


###############################################################################
######### projection to y-z plane
ax2 = fig.add_subplot(223)


### loop over particles
for nnPart in range(nPartPlot): # range(nPart):
	ax2.scatter(pos_vs_t[:,1,nnPart],pos_vs_t[:,2,nnPart],\
		marker='.', s=1) 
	## green start point of trajectory
	ax2.scatter(pos_vs_t[0,1,nnPart],pos_vs_t[0,2,nnPart],\
		marker='o', s=10, color='C2') 	
	## red end point of trajectory
	ax2.scatter(pos_vs_t[nTplot-1,1,nnPart],pos_vs_t[nTplot-1,2,nnPart],\
		marker='o', s=10, color='C3') 
	
ax2.set_xlabel('$y$-position (Px)')
ax2.set_ylabel('$z$-position (Px)')
ax2.axis('equal')
ax2.set_aspect('equal','box')
ax2.set_xticks(np.linspace(0,500,3))
ax2.set_yticks(np.linspace(0,1000,5))
ax2.set_xlim([0,512])
ax2.set_ylim([0,1024])
# plt.show()





###############################################################################
######### projection to x-y plane
ax1 = fig.add_subplot(224)

### loop over particles
for nnPart in range(nPartPlot): # range(nPart):
	ax1.scatter(pos_vs_t[:,0,nnPart],pos_vs_t[:,1,nnPart],\
		marker='.', s=1) 
	## green start point of trajectory
	ax1.scatter(pos_vs_t[0,0,nnPart],pos_vs_t[0,1,nnPart],\
		marker='o', s=10, color='C2') 	
	## red end point of trajectory
	ax1.scatter(pos_vs_t[nTplot-1,0,nnPart],pos_vs_t[nTplot-1,1,nnPart],\
		marker='o', s=10, color='C3') 
	
ax1.set_xlabel('$x$-position (Px)')
ax1.set_ylabel('$y$-position (Px)')
ax1.axis('equal')
ax1.set_aspect('equal','box')
ax1.set_xticks(np.linspace(0,500,6))
ax1.set_yticks(np.linspace(0,500,6))
ax1.set_xlim([0,512])
ax1.set_ylim([0,512])
# plt.show()



#### save plot
fig.savefig(f_meas + dir_out + 'Multipl_trajectory_plot_' +
	'nPart' + str('{:d}'.format(nPartPlot)) + '_' + 
	'nTime' + str('{:d}'.format(nTplot)) + '_' +
	'dynamics_' + str('{:d}'.format(dyN)) +
	'_speed_' + str('{:d}'.format(partSpeedFac)) + '.eps',
	format='eps', bbox_inches='tight') 


plt.close()

import pdb; pdb.set_trace()
