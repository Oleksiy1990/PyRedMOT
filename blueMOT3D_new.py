import numpy as np 
import scipy as sp 
#import matplotlib.pyplot as plt
#import profile
import tables #this is PyTables to use HDF5 for Python 
import sys
import itertools

from scipy.integrate import odeint 

from pyConstants import *
from SrConstants import * 

"""
General equations for the MOT force in 3D. This assumes the approximation of perfect polarization of
the beams
"""

_VIEWPORT_RAD = 19e-3 #radius of the windows in meters (19mm, 38mm diam). This explicitly takes into account the fact that the beam can be clipped by the viewport
_INIT_RAD = 0.1 #This is the radius of the chamber, this doesn't change

def MOT_force(k_laser,linewidth,sat_param,velocity,B_gradient,position,detunings_list):
    """
    This gives the MOT force assuming that polarizations are ideal
    """
    delta_min = [detuning - k_laser*velocity - (muBohr*B_gradient*position/hbar) for detuning in detunings_list]
    delta_plus = [detuning + k_laser*velocity + (muBohr*B_gradient*position/hbar) for detuning in detunings_list]

    #Detunings_list is one with the view of different frequencies in a MOT in mind
    #Note that here by default laser power is not divided up for all the assumed lines. Whatever saturation param is given 
    # will be ised for every line

    force_list = [(hbar*k_laser*linewidth*sat_param/2)*(1/(1+sat_param*(2*delta_min[q]/linewidth)**2) - \
        1/(1+sat_param*(2*delta_plus[q]/linewidth)**2)) for q in range(len(delta_min))] 
    return np.sum(force_list)



def diffeqs_blue(variables,t,params): #We disregard gravity
    (x,vx,y,vy,z,vz) = variables
    (laserPowerX,laserPowerY,laserPowerZ,beamWaistRadX,beamWaistRadY,beamWaistRadZ,B_gradient,detunings_list) = params


    # This assumes circular beams!
    gaussian_Xbeam = (2*laserPowerX/(np.pi*beamWaistRadX**2))*np.exp(-2*(y**2)/beamWaistRadX**2)*np.exp(-2*(z**2)/beamWaistRadX**2) if (y**2+z**2<=_VIEWPORT_RAD**2) else 0
    gaussian_Ybeam = (2*laserPowerY/(np.pi*beamWaistRadY**2))*np.exp(-2*(x**2)/beamWaistRadY**2)*np.exp(-2*(z**2)/beamWaistRadY**2) if (x**2+z**2<=_VIEWPORT_RAD**2) else 0
    gaussian_Zbeam = (2*laserPowerZ/(np.pi*beamWaistRadZ**2))*np.exp(-2*(x**2)/beamWaistRadZ**2)*np.exp(-2*(y**2)/beamWaistRadZ**2) if (x**2+y**2<=_VIEWPORT_RAD**2) else 0


    derivs = [vx,(1/mSr88)*MOT_force(kVecBlue,blueGamma,gaussian_Xbeam/blueIsat,vx,B_gradient,x,detunings_list),\
            vy,(1/mSr88)*MOT_force(kVecBlue,blueGamma,gaussian_Ybeam/blueIsat,vy,B_gradient,y,detunings_list),\
            vz,(1/mSr88)*MOT_force(kVecBlue,blueGamma,gaussian_Zbeam/blueIsat,vz,-2*B_gradient,z,detunings_list)]
            #One needs to be careful with the sign in front of 2B_gradient in the z_dir

    # Quadrupole field from: 
    # https://www2.physics.ox.ac.uk/sites/default/files/2013-01-19/minsung_pdf_16672.pdf eq. 2.40
    return derivs



# Simulation parameters



# num_redCombLines = 50
# detunings_red = np.linspace(-2*np.pi*200*10**3,-2*np.pi*5*10**6,num_redCombLines)


#Put only integer number of milliwatts, not fractions
#bluePowerX = bluePowerY = 15*10**-3 
#bluePowerZ = 2*10**-3
#Put only integer number of millimeters, not fractions
#blueRadX = blueRadY = 10*10**-3
#blueRadZ = 10*10**-3





tStop = 0.3
t = np.linspace(0., tStop, 10**5)

#Parameters that are not changed in the simulation





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




def simulation_function(initial_conditions,beam_radii,powX,powY,powZ,grad,detunings):
    """
    only for all three beams of equal radii
    """
    for radX in beam_radii:
        radY = radZ = radX
        parameters_blue = [powX,powY,powZ,radX,radY,radZ,grad,detunings]

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
            tbl_data.flush()

        descr_save["detuning"] = detunings[0]
        descr_save["gradient"] = grad
        descr_save["init_position_away"] = _INIT_RAD
        descr_save["PowerX"] = bluePowerX
        descr_save["PowerY"] = bluePowerY
        descr_save["PowerZ"] = bluePowerZ
        descr_save["WaistRadX"] = radX
        descr_save["WaistRadY"] = radY
        descr_save["WaistRadZ"] = radZ
        descr_save.append()
        tbl_descriptions.flush()

        file_save.close()

if __name__ == '__main__':
    
    #Parameters that are not changed in the simulation
    init_angle = 30*(np.pi/180) 
    detunings_blue = [-blueGamma]
    blueGradientGcm = 24 # because in the experiment we get now 48 G/cm in z-dir, meaning 24 G/cm horiz.
    blueGradient = 0.01*blueGradientGcm # T/m 

    #NOTE: diam of the last Zeeman slower tube is 21 mm 

    #x_shifts = 0#np.linspace(0,13.86e-3,5)
    #y_shifts = 0#np.linspace(0,24e-3,5)

    #NOTE: in this approach to simulation, I only assume that the atoms start at the center of 
    # the Zeeman slower tube, but have different velocities, not only flying straight 
    # I don't simulate in this code the effects of starting not at the center of the short ZS tube
    initx = -_INIT_RAD*np.sin(init_angle)
    inity = -_INIT_RAD*np.cos(init_angle)
    initz = 0
    inits_pos = np.array([initx,inity,initz])

    init_speed_xy = 30 #m/s
    vel_angle = 30 #deg, this is the angle in the coord sys with respect to the beams
    initvx = init_speed_xy*np.sin(vel_angle*np.pi/180) + np.linspace(-7,7,7)
    initvy = init_speed_xy*np.cos(vel_angle*np.pi/180) + np.linspace(-7,7,7)
    initvz = np.linspace(-5,5,7) 
    inits_vel = np.array(list(itertools.product(initvx,initvy,initvz)))


    inits = np.array([[inits_pos[0],v[0],inits_pos[1],v[1],inits_pos[2],v[2]] for v in inits_vel])
    
    #sys.exit(0)

    
    radiiXYZ = np.array([8,11,14,17])*1e-3 

    #Put only integer number of milliwatts, not fractions
    bluePowerX = bluePowerY = 10*10**-3 
    bluePowerZ = 2*10**-3
    simulation_function(inits,radiiXYZ,bluePowerX,bluePowerY,bluePowerZ,blueGradient,detunings_blue)


    bluePowerX = bluePowerY = 8*10**-3 
    bluePowerZ = 5*10**-3
    simulation_function(inits,radiiXYZ,bluePowerX,bluePowerY,bluePowerZ,blueGradient,detunings_blue)


    bluePowerX = bluePowerY = 7*10**-3 
    bluePowerZ = 7*10**-3
    simulation_function(inits,radiiXYZ,bluePowerX,bluePowerY,bluePowerZ,blueGradient,detunings_blue)




