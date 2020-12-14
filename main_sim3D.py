###################################################################################################
# 3D simulation of particles.  
#
####
# Input:
# nPart :: number of particles
# nTime :: number of time steps
# dyTypes :: different types of dynamcis (directed motion, constant speed, Brownian)
# dimX, dimY, dimZ :: Dimensions of simulation box
#
# Ouput:
# 3d datafiles of particle positions for each timestep.
###################################################################################################


import numpy as np
import time

import module_sim3D # contains functions to create fake images
# import module_velocity_distr_3D # contains functions to create different velocity distributions


#### start run time measurement
start_time = time.time()


##################################################################################################
### set system parameters
depthRatio = 2.
BoxX = 532 ### 512 px + frame for particles of diameter 25px apparing step by step with periodic boundary
BoxY = 532
BoxZ = np.uint32(512*depthRatio)
BoxDim = [BoxX,BoxY,BoxZ] # dimensions of simulation box
nTime = 20 # 000 # number of time steps
dt = 1.0 # step size

#### type of particle dynamics #
####(0 = directed motion, 1 = const speed rand directions, 2 = Max Boltz, 3 = box distr)
#### 4 = sedimentation == directed just downwards
#### 5 = sedimentation with velocity distribution on particle ensemble
#### 6 = sedimentation vel. distr on particle ensemble + overlayed brownian motion
dyN = 0 
nPart = 20 # number of particles (113445 corresponds to packing fraction of 0.4)
speed = 2 ## adjust mean speed of particles between timesteps, displacement in Px
pixel_pitch = 0.009 ## mm/px (scaling factor from mm to px)

#################################################################################
############# output folders for 
############# nPart :: number of particles
############# dyN :: type of dynamics	
	
#### folder of this run
dir_add = ('Sim_output/')




### create folders for output
dir_pos_data, dir_posCon_data, \
dir_vel_data, dir_dynamics = \
module_sim3D.func_dir_out(dir_add, nPart, nTime, dyN, speed)

##################################################
## save simulation parameters
module_sim3D.save_sim_params(dir_add, BoxDim, nTime, dt, dyN, nPart, speed)



### particle dynamics
module_sim3D.particle_dynamics\
	(BoxDim, nPart, nTime, dyN, speed, dt,\
	dir_pos_data, dir_posCon_data, dir_vel_data, dir_dynamics,\
	start_time)




###################################################################################################
##########  END OF PROGRAM 
total_time = time.time() - start_time
tot_min, tot_sec = divmod(total_time, 60)
tot_hrs, tot_min = divmod(tot_min,60)
print ('total time: %d:%02d:%02d' % (tot_hrs,tot_min,tot_sec))


