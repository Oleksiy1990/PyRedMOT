import numpy as np 
import scipy as sp 
import matplotlib.pyplot as plt
#import profile
import random as r

#import tables #this is PyTables to use HDF5 for Python 
import sys
import os

from scipy.integrate import odeint, quad
#from numpy.random import random_sample, uniform
#from mpl_toolkits.mplot3d import Axes3D

from pyConstants import *
from SrConstants import * 

cos = np.cos
sin = np.sin
pi = np.pi

SIDE_VIEWPORT_RAD = 19e-3 #19 mm
VERTICAL_VIEWPORT_RAD = 24e-3


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


    gaussian_Xbeam = np.exp((-2*(y**2+z**2))/beamWaistRadX**2)*(2*laserPowerX)/(np.pi*(beamWaistRadX**2)) if ((y**2+z**2) <= SIDE_VIEWPORT_RAD**2) else 0
    gaussian_Ybeam = np.exp((-2*(x**2+z**2))/beamWaistRadY**2)*(2*laserPowerY)/(np.pi*(beamWaistRadY**2)) if ((x**2+z**2) <= SIDE_VIEWPORT_RAD**2) else 0
    gaussian_Zbeam = np.exp((-2*(x**2+y**2))/beamWaistRadZ**2)*(2*laserPowerZ)/(np.pi*(beamWaistRadZ**2)) if ((x**2+y**2) <= VERTICAL_VIEWPORT_RAD**2) else 0

    # No clue what it was supposed to do
    #phi = uniform(0,np.pi)
    #theta = uniform(0,2*np.pi)


    derivs = [vx,(1/mSr88)*MOT_force(k_vector,Gamma,gaussian_Xbeam/Isat,vx,B_gradient,x,detuning),\
            vy,(1/mSr88)*MOT_force(k_vector,Gamma,gaussian_Ybeam/Isat,vy,B_gradient,y,detuning),\
            vz,(1/mSr88)*MOT_force(k_vector,Gamma,gaussian_Zbeam/Isat,vz,2*B_gradient,z,detuning) - g_Earth]

    # Quadrupole field from: 
    # https://www2.physics.ox.ac.uk/sites/default/files/2013-01-19/minsung_pdf_16672.pdf eq. 2.40
    return derivs

# The functions above define the necessary and general ingredients for simulating MOT dynamics

#================================================

# Simulation parameters

sim_type_counter = 0

while True:

    if sim_type_counter == 0:

        detuning_blue_gammafraction = 0.8
        detuning_blue = -blueGamma*detuning_blue_gammafraction
        bluePowerX = bluePowerY = 9.5e-3 
        bluePowerZ = 2.0e-3
        blueRadX = blueRadY = 15e-3
        blueRadZ = 1e-3
        blueGradient = -0.4 # T/m = 45 G/cm NOTE! Due to sign conventions, this must be written as negative, otherwise equations fail

    elif sim_type_counter == 1:

        detuning_blue_gammafraction = 0.8
        detuning_blue = -blueGamma*detuning_blue_gammafraction
        bluePowerX = bluePowerY = 9.25e-3 
        bluePowerZ = 2.5e-3
        blueRadX = blueRadY = 15e-3
        blueRadZ = 3e-3
        blueGradient = -0.4 

    elif sim_type_counter == 2:

        detuning_blue_gammafraction = 0.8
        detuning_blue = -blueGamma*detuning_blue_gammafraction
        bluePowerX = bluePowerY = 9.0e-3 
        bluePowerZ = 3e-3
        blueRadX = blueRadY = 15e-3
        blueRadZ = 5e-3
        blueGradient = -0.4 

    elif sim_type_counter == 3:

        detuning_blue_gammafraction = 0.8
        detuning_blue = -blueGamma*detuning_blue_gammafraction
        bluePowerX = bluePowerY = 8.75e-3 
        bluePowerZ = 3.5e-3
        blueRadX = blueRadY = 15e-3
        blueRadZ = 7e-3
        blueGradient = -0.4 

    elif sim_type_counter == 4:

        detuning_blue_gammafraction = 0.8
        detuning_blue = -blueGamma*detuning_blue_gammafraction
        bluePowerX = bluePowerY = 8.5e-3 
        bluePowerZ = 4.0e-3
        blueRadX = blueRadY = 15e-3
        blueRadZ = 9e-3
        blueGradient = -0.4 

    elif sim_type_counter == 5:

        detuning_blue_gammafraction = 0.8
        detuning_blue = -blueGamma*detuning_blue_gammafraction
        bluePowerX = bluePowerY = 8.25e-3 
        bluePowerZ = 4.5e-3
        blueRadX = blueRadY = 15e-3
        blueRadZ = 11e-3
        blueGradient = -0.4 

    elif sim_type_counter == 6:

        detuning_blue_gammafraction = 0.8
        detuning_blue = -blueGamma*detuning_blue_gammafraction
        bluePowerX = bluePowerY = 8.0e-3 
        bluePowerZ = 5.0e-3
        blueRadX = blueRadY = 15e-3
        blueRadZ = 13e-3
        blueGradient = -0.4 

    elif sim_type_counter == 7:

        detuning_blue_gammafraction = 0.8
        detuning_blue = -blueGamma*detuning_blue_gammafraction
        bluePowerX = bluePowerY = 7.75e-3 
        bluePowerZ = 5.5e-3
        blueRadX = blueRadY = 15e-3
        blueRadZ = 15e-3
        blueGradient = -0.4

    elif sim_type_counter == 8:

        detuning_blue_gammafraction = 0.8
        detuning_blue = -blueGamma*detuning_blue_gammafraction
        bluePowerX = bluePowerY = 7.5e-3 
        bluePowerZ = 6.0e-3
        blueRadX = blueRadY = 15e-3
        blueRadZ = 17e-3
        blueGradient = -0.4

    elif sim_type_counter == 9:

        detuning_blue_gammafraction = 0.8
        detuning_blue = -blueGamma*detuning_blue_gammafraction
        bluePowerX = bluePowerY = 7.25e-3 
        bluePowerZ = 6.5e-3
        blueRadX = blueRadY = 15e-3
        blueRadZ = 19e-3
        blueGradient = -0.4

    elif sim_type_counter == 10:

        detuning_blue_gammafraction = 0.8
        detuning_blue = -blueGamma*detuning_blue_gammafraction
        bluePowerX = bluePowerY = 7.0e-3 
        bluePowerZ = 7.0e-3
        blueRadX = blueRadY = 15e-3
        blueRadZ = 21e-3
        blueGradient = -0.4

    elif sim_type_counter == 11:

        detuning_blue_gammafraction = 0.8
        detuning_blue = -blueGamma*detuning_blue_gammafraction
        bluePowerX = bluePowerY = 9.0e-3 
        bluePowerZ = 3.0e-3
        blueRadX = blueRadY = 15e-3
        blueRadZ = 1e-3
        blueGradient = -0.4

    elif sim_type_counter == 12:

        detuning_blue_gammafraction = 0.8
        detuning_blue = -blueGamma*detuning_blue_gammafraction
        bluePowerX = bluePowerY = 9.0e-3 
        bluePowerZ = 3.0e-3
        blueRadX = blueRadY = 15e-3
        blueRadZ = 3e-3
        blueGradient = -0.4

    elif sim_type_counter == 13:

        detuning_blue_gammafraction = 0.8
        detuning_blue = -blueGamma*detuning_blue_gammafraction
        bluePowerX = bluePowerY = 9.0e-3 
        bluePowerZ = 3.0e-3
        blueRadX = blueRadY = 15e-3
        blueRadZ = 5e-3
        blueGradient = -0.4

    elif sim_type_counter == 14:

        detuning_blue_gammafraction = 0.8
        detuning_blue = -blueGamma*detuning_blue_gammafraction
        bluePowerX = bluePowerY = 9.0e-3 
        bluePowerZ = 3.0e-3
        blueRadX = blueRadY = 15e-3
        blueRadZ = 7e-3
        blueGradient = -0.4

    elif sim_type_counter == 15:

        detuning_blue_gammafraction = 0.8
        detuning_blue = -blueGamma*detuning_blue_gammafraction
        bluePowerX = bluePowerY = 9.0e-3 
        bluePowerZ = 3.0e-3
        blueRadX = blueRadY = 15e-3
        blueRadZ = 9e-3
        blueGradient = -0.4

    elif sim_type_counter == 16:

        detuning_blue_gammafraction = 0.8
        detuning_blue = -blueGamma*detuning_blue_gammafraction
        bluePowerX = bluePowerY = 9.0e-3 
        bluePowerZ = 3.0e-3
        blueRadX = blueRadY = 15e-3
        blueRadZ = 11e-3
        blueGradient = -0.4

    elif sim_type_counter == 17:

        detuning_blue_gammafraction = 0.8
        detuning_blue = -blueGamma*detuning_blue_gammafraction
        bluePowerX = bluePowerY = 9.0e-3 
        bluePowerZ = 3.0e-3
        blueRadX = blueRadY = 15e-3
        blueRadZ = 13e-3
        blueGradient = -0.4

    elif sim_type_counter == 18:

        detuning_blue_gammafraction = 0.8
        detuning_blue = -blueGamma*detuning_blue_gammafraction
        bluePowerX = bluePowerY = 9.0e-3 
        bluePowerZ = 3.0e-3
        blueRadX = blueRadY = 15e-3
        blueRadZ = 15e-3
        blueGradient = -0.4

    elif sim_type_counter == 19:

        detuning_blue_gammafraction = 0.8
        detuning_blue = -blueGamma*detuning_blue_gammafraction
        bluePowerX = bluePowerY = 9.0e-3 
        bluePowerZ = 3.0e-3
        blueRadX = blueRadY = 15e-3
        blueRadZ = 17e-3
        blueGradient = -0.4

    elif sim_type_counter == 20:

        detuning_blue_gammafraction = 0.8
        detuning_blue = -blueGamma*detuning_blue_gammafraction
        bluePowerX = bluePowerY = 9.0e-3 
        bluePowerZ = 3.0e-3
        blueRadX = blueRadY = 15e-3
        blueRadZ = 19e-3
        blueGradient = -0.4

    elif sim_type_counter > 20:
        sys.exit(0)

    captured_atoms = []
    counter = 0

    simulation_root_dir = "/Users/oleksiy/Desktop/PythonCode/PyRedMOT/resultsBlueMOT_new"
    simulation_directory_name = "sim%i"%sim_type_counter
    simulation_dir_fullpath = simulation_root_dir+"/"+simulation_directory_name



    #This is to prevent simulations from being overwritten
    try:
        os.mkdir(simulation_dir_fullpath)
    except:
        print("Cannot make directory %s for simulaton. Check the name, it probably already exists"%simulation_dir_fullpath)
        sys.exit(0) 

    initial_speeds = np.linspace(5,150,50)

    for init_speed in initial_speeds:

        print("In simulation type %i"%sim_type_counter)

        print("Loop %i out of %i"%(counter,initial_speeds.size))
        counter += 1

        parameters_blue = [kVecBlue,blueGamma,blueIsat,bluePowerX,bluePowerY,bluePowerZ, \
                            blueRadX,blueRadY,blueRadZ,blueGradient,detuning_blue]


        tStop = 0.05
        t = np.linspace(0., tStop, 2e5)


        num_init_angles = 100
        low_end_theta = 90 #Angles in degrees for clarity here
        high_end_theta = 95
        low_end_phi = 25
        high_end_phi = 35

        # the atom starts in the center of the short ZS tube, right at its exit
        inits_blue = [(-0.0866,init_speed*sin(r.uniform(low_end_theta,high_end_theta)*pi/180)*cos(r.uniform(low_end_phi,high_end_phi)*pi/180),\
            -0.05,init_speed*sin(r.uniform(low_end_theta,high_end_theta)*pi/180)*sin(r.uniform(low_end_phi,high_end_phi)*pi/180),\
            0,init_speed*cos(r.uniform(low_end_theta,high_end_theta)*pi/180)) for q in range(num_init_angles)]

        solutions_blue_MOT = [odeint(diffeqs_MOT, inits, t, args=(parameters_blue,),mxstep=10**8) for inits in inits_blue]

        final_positions = np.array([np.sqrt(z[-1,0]**2 + z[-1,2]**2 + z[-1,4]**2) for z in solutions_blue_MOT])
        #final_velocities = np.array([np.sqrt(v[-1,1]**2 + v[-1,3]**2 + v[-1,5]**2) for v in solutions_blue_MOT])

        num_captured = final_positions[np.where(final_positions < 5e-3)].size
        captured_atoms.append(num_captured)

    #turn them into a numpy array for easier saving and calculate the fraction
    captured_atoms_fraction = np.array(captured_atoms)/num_init_angles

    file_description = open(simulation_dir_fullpath+"/"+"simInfo.txt","w")
    #file_results = open(simulation_dir_fullpath+"/"+"results.txt","w")

    file_description.write("Data going into the simulation in this folder: \n")
    file_description.write("Detuning: -{:2f} blueGamma \n power X-beam: {:.4f} W \n power Y-beam: {:.4f} W \n power Z-beam: {:.4f} W \n".format(detuning_blue_gammafraction,bluePowerX,bluePowerY,bluePowerZ))
    file_description.write("Beam radius W_x: {:.4f} m \n Beam radius W_y: {:.4f} m \n Beam radius W_z: {:.4f} m \n".format(blueRadX,blueRadY,blueRadZ))
    file_description.write("B-field grad in radial directions: {:.4f} T/m \n".format(blueGradient))
    file_description.write("Initial position: (x,y,z) = (-0.0866,-0.05,0) m , which is the center at the end of the short ZS tube \n")
    file_description.write("Initial velocities: (vx,vy,vz) = (init_speed*sin(theta)cos(phi),init_speed*sin(theta)sin(phi),init_speed*cos(theta)) m/s \n"+\
                             "where theta is taken from a uniform random distribution [{:d} deg, {:d} deg] and phi from a uniform random distribution [{:d} deg, {:d} deg] \n".format(low_end_theta,high_end_theta,low_end_phi,high_end_phi))

    file_description.close()

    np.savetxt(simulation_dir_fullpath+"/"+"results.txt", (initial_speeds,captured_atoms_fraction),header="Row 0: initial speeds of the atoms [m/s], Row 1: fraction of atoms captured (between 0 and 1), out of {:d}".format(num_init_angles))

    plt.figure(1)
    plt.subplot(1,1,1)
    plt.scatter(initial_speeds,captured_atoms_fraction)
    #plt.title("Final positions")
    plt.xlabel("Initial speed [m/s]")
    plt.ylabel("Fraction of captured atoms")
    plt.savefig(simulation_dir_fullpath+"/"+"figure%i.png"%sim_type_counter)
    plt.savefig("/Users/oleksiy/Desktop/PythonCode/PyRedMOT/resultsBlueMOT_new/Figures/figure%i.png"%sim_type_counter)
    #plt.show()

    sim_type_counter += 1



#print(np.array(final_positions) - np.array(final_velocities))

# plt.scatter(detuning_fraction,captured_atoms)
# plt.show()
# sys.exit(0)

# plt.figure(1)
# plt.subplot(2,1,1)
# plt.scatter(initial_speeds,final_positions)
# #plt.title("Final positions")
# plt.xlabel("Initial speed [m/s]")
# plt.ylabel("Final position [m]")
# plt.ylim(0,2)
# plt.subplot(2,1,2)
# plt.scatter(initial_speeds,final_velocities)
# #plt.title("Final speeds")
# plt.xlabel("Initial speed [m/s]")
# plt.ylabel("Final speed [m/s]")
# plt.ylim(0,2)

# plt.show()



### =============== Tests for a sinle initial condition, plotting entire trajectories ================
# parameters_blue = [kVecBlue,blueGamma,blueIsat,bluePowerX,bluePowerY,bluePowerZ, \
#                         blueRadX,blueRadY,blueRadZ,blueGradient,detuning_blue]


# tStop = 0.05
# t = np.linspace(0., tStop, 2e5)

# speed = 100
# inits_blue_test = (-1*0.0866,speed*sin(91*pi/180)*cos(31*pi/180),-1*0.05,speed*sin(91*pi/180)*sin(31*pi/180),0,speed*cos(91*pi/180))

# solution_blue = odeint(diffeqs_MOT, inits_blue_test, t, args=(parameters_blue,),mxstep=10**8)
# print(solution_blue)
    
# # plt.plot(t,np.sqrt(solution_blue[:,0]**2 + solution_blue[:,2]**2 + solution_blue[:,4]**2),label="pos")
# # plt.plot(t,np.sqrt(solution_blue[:,1]**2 + solution_blue[:,3]**2 + solution_blue[:,5]**2))
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



# 


# fig = plt.figure()
# ax = fig.gca(projection='3d')
# x = solution_blue[:,0]
# y = solution_blue[:,2]
# z = solution_blue[:,4]
# ax.plot(x, y, z, label='parametric curve')
# ax.legend()

# plt.show()
# sys.exit(0)



