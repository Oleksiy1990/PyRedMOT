import numpy as np 
import scipy as sp 
import matplotlib.pyplot as plt
import profile
import h5py
import tables #this is PyTables to use HDF5 for Python 
import sys
import itertools

from scipy.integrate import odeint 
from numpy.random import random_sample, uniform
from mpl_toolkits.mplot3d import Axes3D

from pyConstants import *
from SrConstants import * 

"""
General equations for the MOT force in 3D. This assumes the approximation of perfect polarization of
the beams
"""

viewport_rad = 19e-3 #radius of the windows in meters (19mm, 38mm diam). 
					#This explicitly takes into account the fact that the beam can be clipped by the viewport

def MOT_force(k_laser,linewidth,sat_param,velocity,B_gradient,position,detunings_list):
    delta_min = [detuning - k_laser*velocity - (muBohr*B_gradient*position/hbar) for detuning in detunings_list]
    delta_plus = [detuning + k_laser*velocity + (muBohr*B_gradient*position/hbar) for detuning in detunings_list]
    force_list = [(hbar*k_laser*linewidth*sat_param/2)*(1/(1+sat_param*(2*delta_min[q]/linewidth)**2) - \
        1/(1+sat_param*(2*delta_plus[q]/linewidth)**2)) for q in range(len(delta_min))]
    return np.sum(force_list)



def diffeqs_blue(variables,t,params): #We disregard gravity
    (x,vx,y,vy,z,vz) = variables
    (laserPowerX,laserPowerY,laserPowerZ,beamWaistRadX,beamWaistRadY,beamWaistRadZ,B_gradient,detunings_list) = params


    # This assumes circular beams!
    gaussian_Xbeam = (2*laserPowerX/(np.pi*beamWaistRadX**2))*np.exp(-2*(y**2)/beamWaistRadX**2)*np.exp(-2*(z**2)/beamWaistRadX**2) if (y**2+z**2<=viewport_rad**2) else 0
    gaussian_Ybeam = (2*laserPowerY/(np.pi*beamWaistRadY**2))*np.exp(-2*(x**2)/beamWaistRadY**2)*np.exp(-2*(z**2)/beamWaistRadY**2) if (x**2+z**2<=viewport_rad**2) else 0
    gaussian_Zbeam = (2*laserPowerZ/(np.pi*beamWaistRadZ**2))*np.exp(-2*(x**2)/beamWaistRadZ**2)*np.exp(-2*(y**2)/beamWaistRadZ**2) if (x**2+y**2<=viewport_rad**2) else 0


    derivs = [vx,(1/mSr88)*MOT_force(kVecBlue,blueGamma,gaussian_Xbeam/blueIsat,vx,B_gradient,x,detunings_list),\
            vy,(1/mSr88)*MOT_force(kVecBlue,blueGamma,gaussian_Ybeam/blueIsat,vy,B_gradient,y,detunings_list),\
            vz,0] 

            # we set the force in z-dir to 0 to make it 2D and see how it behaves in that direction 
            

    # Quadrupole field from: 
    # https://www2.physics.ox.ac.uk/sites/default/files/2013-01-19/minsung_pdf_16672.pdf eq. 2.40
    return derivs



# Simulation parameters

detunings_blue = [-blueGamma]


#Put only integer number of milliwatts, not fractions
bluePowerX = bluePowerY = 15*10**-3 
bluePowerZ = 2*10**-3
#Put only integer number of millimeters, not fractions
# blueRadX = blueRadY = 10*10**-3
# blueRadZ = 10*10**-3
blueGradientGcm = 20 #put an integer number here
blueGradient = 0.01*blueGradientGcm # T/m 


# parameters_blue = [bluePowerX,bluePowerY,bluePowerZ,blueRadX,blueRadY,blueRadZ,blueGradient,detunings_blue]



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


class Results(tables.IsDescription):
    x_init = tables.Float64Col()   
    vx_init = tables.Float64Col()
    y_init = tables.Float64Col()
    vy_init = tables.Float64Col()
    z_init = tables.Float64Col()
    vz_init = tables.Float64Col()
    x_final = tables.Float64Col()   
    vx_final = tables.Float64Col()
    y_final = tables.Float64Col()
    vy_final = tables.Float64Col()
    z_final = tables.Float64Col()
    vz_final = tables.Float64Col()
   

class Descriptions(tables.IsDescription):
    detuning = tables.Float64Col()
    gradient = tables.Float64Col()
    init_position_away = tables.Float64Col()
    PowerX = tables.Float64Col()
    PowerY = tables.Float64Col()
    PowerZ = tables.Float64Col() 
    WaistRadX = tables.Float64Col()
    WaistRadY = tables.Float64Col()
    WaistRadZ = tables.Float64Col()



# grp_descriptions = file_save.create_group(grp_simulation,"Descriptions",\
#     title="The input parameters of the simulation are given here, units are SI")

radiiXY = np.array([5,8,11,14,17])*1e-3
for radX in radiiXY:
	radY = radZ = radX
	parameters_blue = [bluePowerX,bluePowerY,bluePowerZ,radX,radY,radZ,blueGradient,detunings_blue]

	file_save = tables.open_file("resultsBlueMOT/pX%.ipY%.ipZ%.igrad%.i.hdf5"%(int(bluePowerX*1e3),int(bluePowerY*1e3),int(bluePowerZ*1e3),\
		blueGradientGcm),mode="a",title= "Blue MOT simulation, detuning -1*gamma")
	grp_simulation = file_save.create_group("/","radX%.iradY%.iradZ%.i"%(int(radX*1e3),int(radY*1e3),int(radZ*1e3)),\
		title="Different beam radii [mm]")
	tbl_data = file_save.create_table(grp_simulation,"Data",Results)
	tbl_descriptions = file_save.create_table(grp_simulation,"Descriptions",Descriptions)
	data_save = tbl_data.row
	descr_save = tbl_descriptions.row
	
	for num,ic in enumerate(inits):
		print("Solving ",num," for rad ",radX)
		solution_blue = odeint(diffeqs_blue, ic, t, args=(parameters_blue,),mxstep=10**8)
		
		data_save["x_init"] = solution_blue[0,0]
		data_save["vx_init"] = solution_blue[0,1]
		data_save["y_init"] = solution_blue[0,2]
		data_save["vy_init"] = solution_blue[0,3]
		data_save["z_init"] = solution_blue[0,4]
		data_save["vz_init"] = solution_blue[0,5]
		data_save["x_final"] = solution_blue[-1,0]
		data_save["vx_final"] = solution_blue[-1,1]
		data_save["y_final"] = solution_blue[-1,2]
		data_save["vy_final"] = solution_blue[-1,3]
		data_save["z_final"] = solution_blue[-1,4]
		data_save["vz_final"] = solution_blue[-1,5]
		data_save.append()

	descr_save["detuning"] = detunings_blue[0]
	descr_save["gradient"] = blueGradient
	descr_save["init_position_away"] = init_rad
	descr_save["PowerX"] = bluePowerX
	descr_save["PowerY"] = bluePowerY
	descr_save["PowerZ"] = bluePowerZ
	descr_save["WaistRadX"] = radX
	descr_save["WaistRadY"] = radY
	descr_save["WaistRadZ"] = radZ
	descr_save.append()

	tbl_data.flush()
	tbl_descriptions.flush()

	file_save.close()


sys.exit(0)



