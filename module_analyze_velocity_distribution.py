##########################################
## Modules for analysis of velocity distributions in simulations
##########################################


import numpy as np

import os

# from matplotlib import pyplot as plt
import matplotlib as mpl
mpl.use('Agg') ## needed to create plots when running via ssh
import matplotlib.pyplot as plt 

### plot settings - latex look
font = {'family' : 'serif',
 	'weight' : 'normal',
 	'size' : 18}
plt.rc('text', usetex=True)
plt.rc('font', **font) # family='serif')
plt.rc


#### create output folders
## Input: directory of simulation data
def func_dir_out(dir_Sim):
	dir_out = dir_Sim + 'output_velocity_distribution_ver2/'
	try:
		os.stat(dir_out)
	except:
		os.mkdir(dir_out)
	
	return dir_out


##### create plot of velocity distributions for all velocity components
def plot_vel_distr(nPart, startTime, endTime, h_marker_fac, nam_out, vx, vy, vz, v_mean):
	### mean velocity (=0 for random motion)
	vx_mean = np.mean(vx)
	vy_mean = np.mean(vy)
	vz_mean = np.mean(vz)

	
	
		
	############# check which quantity I want to measure and how 
	vtot_mean = np.mean(v_mean) ## = np.mean(np.sqrt(vx**2. + vy**2. + vz**2.))
	# vtot_mean = np.sqrt(vx_mean**2. + vy_mean**2. + vz_mean**2.)
	# vtot_mean = np.sqrt(np.mean(vx**2.) + np.mean(vy_mean**2.) + np.mean(vz_mean**2.))
	# vtot_mean = np.sqrt(np.mean(vx**2.) + np.mean(vy**2.) + np.mean(vz**2.))

	### mean of absolute speed
	vx_abs_mean = np.mean(np.abs(vx))
	vy_abs_mean = np.mean(np.abs(vy))
	vz_abs_mean = np.mean(np.abs(vz))

	vx_std = np.std(vx)
	vy_std = np.std(vy)
	vz_std = np.std(vz) 
	vtot_std = np.std(v_mean)
	
	h_marker = (endTime-startTime) * h_marker_fac	
	
	##########################################################################################
	########## figure ########################################################################
	# import pdb; pdb.set_trace()	

	fig = plt.figure()
	fig.set_size_inches(12,14)
	
	########## x-velocities
	ax0 = fig.add_subplot(221) # the first row of graphs
	ax0.set_title(('$\langle v_x \\rangle$ = ' + str('{:.5f}'.format(vx_mean)) + ', ' + \
		'$\delta v_x$ = ' + str('{:.5f}'.format(vx_std))))
	ax0.set_xlabel('$v_x$ (Px/timestep)')
	ax0.set_ylabel('counts')
	ax0.tick_params(which='both',direction = 'in')
	ax0.xaxis.set_ticks_position('both')
	ax0.yaxis.set_ticks_position('both')
	plt.hist(vx,bins=np.int(np.sqrt(nPart)), color='C0')     
	plt.plot([vx_mean, vx_mean],[0,h_marker], color='C1') 
	plt.plot([vx_mean - vx_std, vx_mean - vx_std], [0,h_marker], color='C2')
	plt.plot([vx_mean + vx_std, vx_mean + vx_std], [0,h_marker], color='C2')

	########## y-velocities
	ax1 = fig.add_subplot(222)
	ax1.set_title(('$\langle v_y \\rangle$ = ' + str('{:.5f}'.format(vy_mean)) + ', ' + \
		'$\delta v_y$ = ' + str('{:.5f}'.format(vy_std))))
	ax1.set_xlabel('$v_y$ (Px/timestep)')
	ax1.set_ylabel('counts')
	ax1.tick_params(which='both',direction = 'in')
	ax1.xaxis.set_ticks_position('both')
	ax1.yaxis.set_ticks_position('both')
	plt.hist(vy,bins=np.int(np.sqrt(nPart)))
	plt.plot([vy_mean, vy_mean],[0,h_marker], color='C1') 
	plt.plot([vy_mean - vy_std, vy_mean - vy_std], [0,h_marker], color='C2')
	plt.plot([vy_mean + vy_std, vy_mean + vy_std], [0,h_marker], color='C2')

	########## z-velocities
	ax2 = fig.add_subplot(223) # the second contains
	ax2.set_title(('$\langle v_z \\rangle$ = ' + str('{:.5f}'.format(vz_mean)) + ', ' + \
		'$\delta v_z$ = ' + str('{:.5f}'.format(vz_std))))
	ax2.set_xlabel('$v_z$ (Px/timestep)')
	ax2.set_ylabel('counts')
	ax2.tick_params(which='both',direction = 'in')
	ax2.xaxis.set_ticks_position('both')
	ax2.yaxis.set_ticks_position('both')
	plt.hist(vz,bins=np.int(np.sqrt(nPart)))     
	plt.plot([vz_mean, vz_mean],[0,h_marker], color='C1') 
	plt.plot([vz_mean - vz_std, vz_mean - vz_std], [0,h_marker], color='C2')
	plt.plot([vz_mean + vz_std, vz_mean + vz_std], [0,h_marker], color='C2')

	########## total-velocites
	ax3 = fig.add_subplot(224)
	ax3.set_title(('$\langle v_{total} \\rangle$ = ' + str('{:.5f}'.format(vtot_mean)) + ', ' + \
		'$\delta v_{tot}$ = ' + str('{:.5f}'.format(vtot_std))))
	ax3.set_xlabel('$v_{total}$ (Px/timestep)')
	ax3.set_ylabel('counts')
	ax3.tick_params(which='both',direction = 'in')
	ax3.xaxis.set_ticks_position('both')
	ax3.yaxis.set_ticks_position('both')
	plt.hist(v_mean,bins=np.int(np.sqrt(nPart)))
	plt.plot([vtot_mean, vtot_mean],[0,h_marker], color='C1') 

	
	########### save figure		
	fig.savefig(nam_out + '.png',
		format='png', bbox_inches='tight')
	plt.close()


	########### save mean and std of velocity distributions
	file_out = open((nam_out + '.dat'), 'w+')
	file_out.write('# <vx>\t <vy>\t <vz>\t <v_tot>\n')
	str_mean_vel = '%.10f %.10f %.10f %.10f\n'%(vx_mean, vy_mean, vz_mean, vtot_mean)
	file_out.write(str_mean_vel)
	file_out.write('\n')
	
	file_out.write('# dvx\t dvy\t dvz\t dv_tot\n')
	str_std_vel = '%.10f %.10f %.10f %.10f\n'%(vx_std, vy_std, vz_std, vtot_std)
	file_out.write(str_std_vel)
	file_out.write('\n')
	
	file_out.write('# <|vx|>\t <|vy|>\t <|vz|>\t <v_tot_abs>\n')
	str_mean_abs_vel = '%.10f %.10f %.10f %.10f\n'%(vx_abs_mean, vy_abs_mean, vz_abs_mean, vtot_mean)
	file_out.write(str_mean_abs_vel)
	
	

