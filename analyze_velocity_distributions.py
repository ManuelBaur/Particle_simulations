######################################
#### velocity distributions from simualtion data
######################################
import numpy as np

import time

import module_analyze_velocity_distribution

#### start run time measurement
start_time = time.time()


### simulation parameters
nPart = 20 ## number of particles
nTime = 20

### name of simultion output folder
Sim_nam = ('Sim_output/')



### directory to simulation output folder
f_Sim = ('' + Sim_nam)
### directory to particle velocity data
f_data = (f_Sim + 'velocites_Part_' + str('{:0>8d}'.format(nPart)) + \
	'_nTime' + str('{:0>4d}'.format(nTime)) + '/')

### create output folder
dir_out = module_analyze_velocity_distribution.func_dir_out(f_Sim)

#### how many timesteps to analyze
lenTime = 10 
startTime = 1 
h_marker_fac = 10 ## factor for height of marker in plot
endTime = startTime + lenTime

### allocate empty vectors for vx, vy, vz, v_mean
vx = []
vy = []
vz = []
v_mean = []

### loop over timesteps
for nnTime in range(startTime, endTime):
	#### time processed
	if nnTime % 10 == 0:
		print ('working on timestep ' + str(startTime + nnTime) +
			' out of ' + str(endTime-startTime))
		elapsed_time = time.time() - start_time
		el_min, el_sec = divmod(elapsed_time, 60)
		el_hrs, el_min = divmod(el_min, 60) 
		print ('elapsed time: %d:%02d:%02d' % (el_hrs,el_min,el_sec))
	
	### read names of velocity files
	nam_data = (f_data + \
		'particle_velocities_nPart' + str('{:0>8d}'.format(nPart)) + \
		'_timestep' + str('{:0>8d}'.format(nnTime)) + '.dat')
	
	##### load all velocity components of all particiels
	# import pdb; pdb.set_trace()
	data_in = np.loadtxt(nam_data)
	vx = np.append(vx, data_in[:,0]) ### x-velocities 
	vy = np.append(vy, data_in[:,1]) ### y-velocities
	vz = np.append(vz, data_in[:,2]) ### z-velocities
	v_mean = np.append(v_mean, data_in[:,3]) ### == sqrt(vx**2 + vy**2 + vz**2) 

##################################################	
##### generate plot of distributions

##### double precision of velocities
nam_out = (dir_out + 'particle_velocity_distribution_nPart' +
		str('{:0>8d}'.format(nPart)) +
		'_nTimestep' + str('{:0>8}'.format(lenTime)))

module_analyze_velocity_distribution.plot_vel_distr\
	(nPart, startTime, endTime, h_marker_fac, nam_out, vx, vy, vz, v_mean)

#### velocites floored to integer values
vx_fl = np.int32(np.floor(vx))
vy_fl = np.int32(np.floor(vy))
vz_fl = np.int32(np.floor(vz))
v_mean_fl = np.sqrt(vx_fl**2 + vy_fl**2 + vz_fl**2)

nam_out_floor = (dir_out + 'particle_velocity_distribution_floored_to_int_nPart' +
		str('{:0>8d}'.format(nPart)) +
		'_nTimestep' + str('{:0>8}'.format(lenTime)))

module_analyze_velocity_distribution.plot_vel_distr\
	(nPart, startTime, endTime, h_marker_fac*30, nam_out_floor, vx_fl, vy_fl, vz_fl, v_mean_fl)

#### velocites rounded to integer values
vx_rd = np.int32(np.round(vx))
vy_rd = np.int32(np.round(vy))
vz_rd = np.int32(np.round(vz))
v_mean_rd = np.sqrt(vx_rd**2 + vy_rd**2 + vz_rd**2)

nam_out_round = (dir_out + 'particle_velocity_distribution_round_to_int_nPart' +
		str('{:0>8d}'.format(nPart)) +
		'_nTimestep' + str('{:0>8}'.format(lenTime)))

module_analyze_velocity_distribution.plot_vel_distr\
	(nPart, startTime, endTime, h_marker_fac*30, nam_out_round, vx_rd, vy_rd, vz_rd, v_mean_rd)


###################################################################################################
##########  END OF PROGRAM 
total_time = time.time() - start_time
tot_min, tot_sec = divmod(total_time, 60)
tot_hrs, tot_min = divmod(tot_min,60)
print ('total time: %d:%02d:%02d' % (tot_hrs,tot_min,tot_sec))


