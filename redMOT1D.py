import numpy as np 
import scipy as sp 
#import matplotlib.pyplot as plt
import tables #this is PyTables to use HDF5 for Python 
#import itertools
import sys
import os

#import multiprocessing as mpr 
from multiprocessing import Process, Lock


from scipy.integrate import odeint 
#from numpy.random import random_sample, uniform
#from mpl_toolkits.mplot3d import Axes3D

from pyConstants import *
from SrConstants import * 


windowRadius = 19e-3 
#this will be the cutoff for the action of the beam, if the atom is outside, the intensity is set to 0 

def MOT_force(k_laser,linewidth,sat_param,velocity,B_gradient,position,detunings_list): 
    """
    detunings must be given as a List data structure, even if there's only a single detuning given
    """ 
    delta_min = [detuning - k_laser*velocity - (muBohr*B_gradient*position/hbar) for detuning in detunings_list]
    delta_plus = [detuning + k_laser*velocity + (muBohr*B_gradient*position/hbar) for detuning in detunings_list]
    force_list = [(hbar*k_laser*linewidth*sat_param/2)*(1/(1+sat_param*(2*delta_min[q]/linewidth)**2) - \
        1/(1+sat_param*(2*delta_plus[q]/linewidth)**2)) for q in range(len(delta_min))]
    return np.sum(force_list)



def diffeqs1D_red(variables,t,params):
    
    (x,vx) = variables
    (power,waistRad,B_gradient,r,gravity,detunings_list) = params
    # give the "gravity" paramter as 1 or 0 to have it ON or OFF

    redIntensity = (2*power/(np.pi*waistRad**2))*np.exp(-2*r**2/waistRad**2)/len(detunings_list) if abs(x) < windowRadius else 0
    #This division is important because it's the whole reason why we need more power in the broadband red MOT

    derivs = [vx,(1/mSr88)*MOT_force(kVecRed,redGamma,redIntensity/redIsat,vx,B_gradient,x,detunings_list) - gravity*9.81] #gravity can be either turned on or off
            # NOTE I removed the minus sign in the z-comp of B-field 
    return derivs

    # Quadrupole field from: 
    # https://www2.physics.ox.ac.uk/sites/default/files/2013-01-19/minsung_pdf_16672.pdf eq. 2.40

# def laguerre_gauss_int(rad,power,waistRad):
#     intensity = (2*power/(np.pi*waistRad**2))*np.exp(-2*rad**2/waistRad**2)
#     return intensity

class Results(tables.IsDescription):
    r_init = tables.Float64Col() 
    x_init = tables.Float64Col()
    vx_init = tables.Float64Col()
    x_final = tables.Float64Col()
    vx_final = tables.Float64Col()
        
class SimParameters(tables.IsDescription):
    redPower = tables.Float64Col()
    redRadius = tables.Float64Col()
    Bgradient = tables.Float64Col()

class Detunings(tables.IsDescription):
    """docstring for Detunings"""
    detunings = tables.Float64Col()


def simulation_func_red(init,params,t,tbl_results,tbl_simparameters,current_init,current_param):

    print("Running init %.i, param %.i"%(current_init,current_param))

    results = tbl_results.row
    parameters = tbl_simparameters.row

    parameters["redPower"] = params[0]
    parameters["redRadius"] = params[1]
    parameters["Bgradient"] = params[2]
    parameters.append()

    tbl_simparameters.flush()

    solution = odeint(diffeqs1D_red, init, t, args=(params,),full_output=False,mxstep=10**4) #it was 10**8 before
    #print(solution[1])
    results["r_init"] = params[3]
    results["x_init"] = solution[0,0]
    results["vx_init"] = solution[0,1]
    results["x_final"] = solution[-1,0]
    results["vx_final"] = solution[-1,1]
    results.append()

    tbl_results.flush()
        

if __name__ == '__main__': 

    num_redCombLines = int(sys.argv[1]) #Give the number of comb lines as a command line argument 
    for redPower in [7e-3,10e-3,13e-3]:
        for redBeamRadius in [2e-3,7e-3,12e-3,17e-3]:

            print("Doing beam rad %.i, power %.i, lines %.i"%(int(redBeamRadius*1e3),int(redPower*1e3),num_redCombLines))

            #num_redCombLines = 250
            detunings_red = np.linspace(-2*np.pi*200*10**3,-2*np.pi*5*10**6,num_redCombLines)
            redGradient = 1.15/100 # T/m
            gravity = 1 #this is given only as 1 or 0 to determine if it's taken into account or not 


            #redPower = 10e-3
            #redBeamRadius = 17e-3 

            params = [[redPower,redBeamRadius,redGradient,r,gravity,detunings_red] for r in np.linspace(0,17e-3,7)]
            inits = [(x,vx) for x in np.linspace(-17e-3,0,10) for vx in np.linspace(-0.9,0.9,10)] #0.9 m/s is the max speed covering >95% atoms at 720 microK
            print("There are %.i params"%len(params))
            print("There are %.i inits"%len(inits))

            tStop = 0.5
            t = np.linspace(0., tStop, 10**4)

            naming_tuple = (int(params[0][0]*1e3), # this gives the power in mW
                int(params[0][1]*1e3), # this gives the beam radius in mm 
                int(params[0][2]*1e4), # this gives the gradient, with the assumption of a comma after the first digit
                len(params[0][5])) # this gives the number of comb lines in broadband MOT


            
            filename_windows = "C:/Users/Oleksiy/Desktop/Code/PyRedMOT/resultsRedMOT/redMOT%.ilinesGravityOn.hdf5"%num_redCombLines
            filename_mac = "/Users/oleksiy/Desktop/PythonCode/PyRedMOT/resultsRedMOT/redMOT%.ilines.hdf5"%num_redCombLines

            if bool(gravity):
                file_save = tables.open_file(filename_windows,mode="a",title= "Red MOT simulation, %.i comb"%num_redCombLines)
            else:
                file_save = tables.open_file(filename_windows,mode="a",title= "Red MOT simulation, %.i comb"%num_redCombLines)
            

            simulation_groupname = "gradient%.iGcm"%int(redGradient*1e4)

            if not file_save.__contains__("/"+simulation_groupname):
                file_save.create_group("/",simulation_groupname)

            table_res_name = "Pow%.imWredRad%.immGrad%.iGcmlines%.i"%naming_tuple
            table_simparams_name = "ParamsPow%.imWredRad%.immGrad%.iGcmlines%.i"%naming_tuple
            array_detunigns_name = "ComblinesPow%.imWredRad%.immGrad%.iGcmlines%.i"%naming_tuple

            if file_save.__contains__("/"+simulation_groupname+"/"+table_res_name):
                file_save.remove_node("/"+simulation_groupname+"/"+table_res_name,recursive=True)
            if file_save.__contains__("/"+simulation_groupname+"/"+table_simparams_name):
                file_save.remove_node("/"+simulation_groupname+"/"+table_simparams_name,recursive=True)
            if file_save.__contains__("/"+simulation_groupname+"/"+array_detunigns_name):
                file_save.remove_node("/"+simulation_groupname+"/"+array_detunigns_name,recursive=True)
            
            tbl_results = file_save.create_table("/gradient%.iGcm"%int(redGradient*1e4),table_res_name,Results,"Results") 
            tbl_simparameters = file_save.create_table("/gradient%.iGcm"%int(redGradient*1e4),table_simparams_name,SimParameters,"Parameters")
            
            arr_detunings = file_save.create_array("/gradient%.iGcm"%int(redGradient*1e4),array_detunigns_name,params[0][5],"Comb lines")
              
            
            [simulation_func_red(init,param,t,tbl_results,tbl_simparameters,ni,np) for ni,init in enumerate(inits) for np,param in enumerate(params)]
            file_save.close()




