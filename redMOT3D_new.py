import numpy as np 
import scipy as sp 
import matplotlib.pyplot as plt
#import profile
import random as r

import tables #this is PyTables to use HDF5 for Python 
import sys

from scipy.integrate import odeint, quad
from numpy.random import random_sample, uniform
from mpl_toolkits.mplot3d import Axes3D

from pyConstants import *
from SrConstants import * 

def MOT_force(k_laser,linewidth,sat_param,velocity,B_gradient,position,detuning):
    #sat_param = I/I_sat
    delta_min = detuning + k_laser*velocity - (muBohr*B_gradient*position/hbar)
    delta_plus = detuning - k_laser*velocity + (muBohr*B_gradient*position/hbar)
    force_max = hbar*k_laser*linewidth/2 # this is assuming infinite intensity
    force_total = force_max*(sat_param/(1+sat_param+(2*delta_plus/linewidth)**2) - \
        sat_param/(1+sat_param+(2*delta_min/linewidth)**2)) 
    #Now the equations for the force look correct 
    return force_total


def diffeqs_MOT(variables,t,params):
    (x,vx,y,vy,z,vz) = variables
    (k_vector,Gamma,Isat,laserPowerX,laserPowerY,laserPowerZ,\
        beamWaistRadX,beamWaistRadY,beamWaistRadZ,B_gradient,detuning) = params


    gaussian_Xbeam = np.exp((-2*(y**2+z**2))/beamWaistRadX**2)*(2*laserPowerX)/(np.pi*(beamWaistRadX**2))
    gaussian_Ybeam = np.exp((-2*(x**2+z**2))/beamWaistRadY**2)*(2*laserPowerY)/(np.pi*(beamWaistRadY**2))
    gaussian_Zbeam = np.exp((-2*(x**2+y**2))/beamWaistRadZ**2)*(2*laserPowerZ)/(np.pi*(beamWaistRadZ**2))

    # No clue what it was supposed to do
    #phi = uniform(0,np.pi)
    #theta = uniform(0,2*np.pi)


    derivs = [vx,(1/mSr88)*MOT_force(k_vector,Gamma,gaussian_Xbeam/Isat,vx,B_gradient,x,detuning),\
            vy,(1/mSr88)*MOT_force(k_vector,Gamma,gaussian_Ybeam/Isat,vy,B_gradient,y,detuning),\
            vz,(1/mSr88)*MOT_force(k_vector,Gamma,gaussian_Zbeam/Isat,vz,2*B_gradient,z,detuning) - g_Earth]

    # Quadrupole field from: 
    # https://www2.physics.ox.ac.uk/sites/default/files/2013-01-19/minsung_pdf_16672.pdf eq. 2.40
    return derivs

# Simulation parameters



# This is the simulation for red, let's disregard the blue for now


detuning_blue = -blueGamma*1
bluePowerX = bluePowerY = 6.5e-3 
bluePowerZ = 6.5e-3
blueRadX = blueRadY = 15e-3
blueRadZ = 15e-3
blueGradient = -0.55 # T/m = 45 G/cm NOTE! Due to sign conventions, this must be written as negative, otherwise equations fail

detuning_fraction = np.linspace(0.5,3.,10)
captured_atoms = []
counter = 0

for fr in detuning_fraction:

	print("Loop %i out of %i"%(counter,detuning_fraction.size))
	counter += 1

	detuning_blue = -blueGamma*fr

	parameters_blue = [kVecBlue,blueGamma,blueIsat,bluePowerX,bluePowerY,bluePowerZ, \
	                    blueRadX,blueRadY,blueRadZ,blueGradient,detuning_blue]


	tStop = 0.1
	t = np.linspace(0., tStop, 1e5)


	num_init_vels = 200
	low_end_vel = 1
	high_end_vel = 30

	inits_blue = [(0,r.uniform(low_end_vel,high_end_vel),0,r.uniform(low_end_vel,high_end_vel),0,-r.uniform(low_end_vel,high_end_vel)) for q in range(num_init_vels)]
	#sys.exit(0)

	#inits_blue_test = (0,5,0,5,0,5)
	solutions_blue_MOT = [odeint(diffeqs_MOT, inits, t, args=(parameters_blue,),mxstep=10**8) for inits in inits_blue]


#solution_blue = odeint(diffeqs_MOT, inits_blue_test, t, args=(parameters_blue,),mxstep=10**8)
#print(solution_blue)

	initial_speeds = np.array([np.sqrt(q[1]**2 + q[3]**2 + q[5]**2) for q in inits_blue])
	final_positions = np.array([np.sqrt(z[-1,0]**2 + z[-1,2]**2 + z[-1,4]**2) for z in solutions_blue_MOT])
	final_velocities = np.array([np.sqrt(v[-1,1]**2 + v[-1,3]**2 + v[-1,5]**2) for v in solutions_blue_MOT])

	num_captured = final_positions[np.where(final_positions < 0.05)].size
	captured_atoms.append(num_captured)
#print(np.array(final_positions) - np.array(final_velocities))

plt.scatter(detuning_fraction,captured_atoms)
plt.show()
sys.exit(0)


plt.figure(1)
plt.subplot(2,1,1)
plt.scatter(initial_speeds,final_positions)
#plt.title("Final positions")
plt.xlabel("Initial speed [m/s]")
plt.ylabel("Final position [m]")
plt.ylim(0,2)
plt.subplot(2,1,2)
plt.scatter(initial_speeds,final_velocities)
#plt.title("Final speeds")
plt.xlabel("Initial speed [m/s]")
plt.ylabel("Final speed [m/s]")
plt.ylim(0,2)

plt.show()


# fig = plt.figure()
# ax = fig.gca(projection='3d')
# x = solution_blue[:,0]
# y = solution_blue[:,2]
# z = solution_blue[:,4]
# ax.plot(x, y, z, label='parametric curve')
# ax.legend()

# plt.show()
# sys.exit(0)

#plt.plot(t,np.sqrt(solution_blue[:,0]**2 + solution_blue[:,2]**2 + solution_blue[:,4]**2),label="pos")
#plt.plot(t,np.sqrt(solution_blue[:,1]**2 + solution_blue[:,3]**2 + solution_blue[:,5]**2))
# plt.plot(t,solution_blue[:,1],label="vx")
# plt.plot(t,solution_blue[:,3],label="vy")
# plt.plot(t,solution_blue[:,5],label="vz")
# plt.legend()
# plt.show()

# plt.plot(t,solution_blue[:,0],label="x")
# plt.plot(t,solution_blue[:,2],label="y")
# plt.plot(t,solution_blue[:,4],label="z")
# plt.legend()
# plt.show()


