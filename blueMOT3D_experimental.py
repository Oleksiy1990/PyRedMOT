import numpy as np 
import scipy as sp 
import matplotlib.pyplot as plt
#import profile
#import h5py
import tables #this is PyTables to use HDF5 for Python 
import sys
import itertools

from scipy.integrate import odeint 
from numpy.random import random_sample, uniform
from mpl_toolkits.mplot3d import Axes3D

from pyConstants import *
from SrConstants import * 


viewport_rad = 19e-3 #radius of the windows in meters (19mm, 38mm diam). This explicitly takes into account the fact that the beam can be clipped by the viewport

def MOT_force(k_laser,linewidth,sat_param,velocity,B_gradient,position,detunings_list):
    delta_min = [detuning - k_laser*velocity - (muBohr*B_gradient*position/hbar) for detuning in detunings_list]
    delta_plus = [detuning + k_laser*velocity + (muBohr*B_gradient*position/hbar) for detuning in detunings_list]
    force_list = [(hbar*k_laser*linewidth*sat_param/2)*(1/(1+sat_param*(2*delta_min[q]/linewidth)**2) - \
        1/(1+sat_param*(2*delta_plus[q]/linewidth)**2)) for q in range(len(delta_min))]
    return np.sum(force_list)



def diffeqs_blue(variables,t,params): #We disregard gravity
    (x,vx) = variables
    (laserPowerX,beamWaistRadX,B_gradient,detunings_list) = params


    # This assumes circular beams!
    gaussian_Xbeam = (2*laserPowerX/(np.pi*beamWaistRadX**2))*np.exp(-2*(y**2)/beamWaistRadX**2)*np.exp(-2*(z**2)/beamWaistRadX**2) if (y**2+z**2<=viewport_rad**2) else 0
    gaussian_Ybeam = (2*laserPowerY/(np.pi*beamWaistRadY**2))*np.exp(-2*(x**2)/beamWaistRadY**2)*np.exp(-2*(z**2)/beamWaistRadY**2) if (x**2+z**2<=viewport_rad**2) else 0
    gaussian_Zbeam = (2*laserPowerZ/(np.pi*beamWaistRadZ**2))*np.exp(-2*(x**2)/beamWaistRadZ**2)*np.exp(-2*(y**2)/beamWaistRadZ**2) if (x**2+y**2<=viewport_rad**2) else 0


    derivs = [vx,(1/mSr88)*MOT_force(kVecBlue,blueGamma,gaussian_Xbeam/blueIsat,vx,B_gradient,x,detunings_list),\
            vy,(1/mSr88)*MOT_force(kVecBlue,blueGamma,gaussian_Ybeam/blueIsat,vy,B_gradient,y,detunings_list),\
            vz,(1/mSr88)*MOT_force(kVecBlue,blueGamma,gaussian_Zbeam/blueIsat,vz,2*B_gradient,z,detunings_list)]
            #One needs to be careful with the sign in front of 2B_gradient in the z_dir

    # Quadrupole field from: 
    # https://www2.physics.ox.ac.uk/sites/default/files/2013-01-19/minsung_pdf_16672.pdf eq. 2.40
    return derivs



# Simulation parameters

detunings_blue = [-blueGamma]
#detunings_red = [-2*np.pi*200*10**3]

# num_redCombLines = 50
# detunings_red = np.linspace(-2*np.pi*200*10**3,-2*np.pi*5*10**6,num_redCombLines)


#Put only integer number of milliwatts, not fractions
bluePowerX = bluePowerY = 15*10**-3 
bluePowerZ = 2*10**-3
#Put only integer number of millimeters, not fractions
blueRadX = blueRadY = 10*10**-3
blueRadZ = 10*10**-3
blueGradientGcm = 20 #put an integer number here
blueGradient = 0.01*blueGradientGcm # T/m = 55 G/cm



tStop = 0.5
t = np.linspace(0., tStop, 10**5)

initz = np.linspace(-12e-3,12e-3,7)
init_rad = 0.1
init_angle = 30*(np.pi/180) #not to change
x_shifts = np.linspace(0,13.86e-3,5)
y_shifts = np.linspace(0,24e-3,5)
initx = -init_rad*np.sin(init_angle)+x_shifts
inity = -init_rad*np.cos(init_angle)+y_shifts
inits_pos_1 = np.array(list(itertools.product(initx,[inity[0]],initz)))
inits_pos_2 = np.array(list(itertools.product([initx[0]],inity[1:],initz)))
inits_pos = np.concatenate((inits_pos_1,inits_pos_2))


init_speed_xy = 40 #m/s
vel_angle = np.linspace(25,35,5) #deg
initvx = init_speed_xy*np.sin(vel_angle*np.pi/180)
initvy = init_speed_xy*np.cos(vel_angle*np.pi/180)
initvz = np.linspace(0,2,4) 
inits_vel = np.array([[init_speed_xy*np.sin(angle*np.pi/180),init_speed_xy*np.cos(angle*np.pi/180),vz] for angle in vel_angle for vz in initvz])


inits = np.array([[p[0],v[0],p[1],v[1],p[2],v[2]] for p in inits_pos for v in inits_vel])
parameters_blue = [bluePowerX,bluePowerY,bluePowerZ,blueRadX,blueRadY,blueRadZ,blueGradient,detunings_blue]

odeint(diffeqs_blue, inits, t, args=(parameters_blue,),mxstep=10**8)

sys.exit(0)

