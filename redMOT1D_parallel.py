import numpy as np 
import scipy as sp 
import matplotlib.pyplot as plt
import tables #this is PyTables to use HDF5 for Python 
import itertools
import sys
import os

import multiprocessing as mpr 



from scipy.integrate import odeint 
from numpy.random import random_sample, uniform
from mpl_toolkits.mplot3d import Axes3D

from pyConstants import *
from SrConstants import * 


#lck = mpr.Lock()

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
    (power,waistRad,B_gradient,r,detunings_list) = params

    redIntensity = (2*power/(np.pi*waistRad**2))*np.exp(-2*r**2/waistRad**2)/len(detunings_list)
    #This division is important because it's the whole reason why we need more power in the broadband red MOT

    derivs = [vx,(1/mSr88)*MOT_force(kVecRed,redGamma,redIntensity/redIsat,vx,B_gradient,x,detunings_list) - 0*9.81]
            # NOTE I removed the minus sign in the z-comp of B-field 
    return derivs

    # Quadrupole field from: 
    # https://www2.physics.ox.ac.uk/sites/default/files/2013-01-19/minsung_pdf_16672.pdf eq. 2.40

# def laguerre_gauss_int(rad,power,waistRad):
#     intensity = (2*power/(np.pi*waistRad**2))*np.exp(-2*rad**2/waistRad**2)
#     return intensity

class Results(tables.IsDescription):
    r_init = tables.Float64Col()   
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


def simulation_func_red(init,params,t,tbl_results,tbl_simparameters):

    print("Getting ready to run the solution loop")


    # print("Process id: ",os.getpid())

    results = tbl_results.row
    parameters = tbl_simparameters.row

    parameters["redPower"] = params[0]
    parameters["redRadius"] = params[1]
    parameters["Bgradient"] = params[2]
    parameters.append()

    tbl_simparameters.flush()

    solution = odeint(diffeqs1D_red, init, t, args=(params,),mxstep=10**8)
    
    results["r_init"] = params[3]
    results["vx_init"] = solution[0,1]
    results["x_final"] = solution[-1,0]
    results["vx_final"] = solution[-1,1]
    results.append()

    tbl_results.flush()


        
        
        
        # #we work with locks so that the file doesn't get corrupted by simultaneous saving
        # lck.acquire()
        # tbl_results.flush()
        # lck.release()
        


        # del tbl_results
        # del solution_red
        # print("Initial conditions %.i out of %.i are done"%(counter,len(initialconds_red)))
        # counter += 1

if __name__ == '__main__': 

    num_redCombLines = 250
    detunings_red = np.linspace(-2*np.pi*200*10**3,-2*np.pi*5*10**6,num_redCombLines)
    redGradient = 1.15/100 # T/m 

    
    redPower = 20e-3
    redBeamRadius = 10e-3 

    params = [[redPower,redBeamRadius,redGradient,r,detunings_red] for r in np.linspace(0,10e-3,10)]
    inits = [(0,vx) for vx in np.linspace(0,1.25,10)]
    [print(init) for init in inits]

    tStop = 0.5
    t = np.linspace(0., tStop, 10**5)

    naming_tuple = (int(params[0][0]*1e3),int(params[0][1]*1e3),int(params[0][2]*1e4),len(params[0][4]))


    

    grp_simulation = file_save.create_group("/","gradient%.iGcm"%int(redGradient*1e4))
    
    tbl_results = file_save.create_table(grp_simulation,"Pow%.imWredRad%.immGrad%.iGcmlines%.i"%naming_tuple,Results,"Results") 
    tbl_simparameters = file_save.create_table(grp_simulation,"ParamsPow%.imWredRad%.immGrad%.iGcmlines%.i"%naming_tuple,SimParameters,"Parameters")
    
    arr_detunings = file_save.create_array(grp_simulation,"ComblinesPow%.imWredRad%.immGrad%.iGcmlines%.i"%naming_tuple,params[0][4],"Comb lines")
      
    
    [simulation_func_red(init,param,t,tbl_results,tbl_simparameters) for init in inits for param in params]
    file_save.close()

    # grp_detunings = file_save.create_group(grp_simulation,"detuningsHz",\
    #     title="Values of freqs of the comb lines, in Hz, red-detuned from resonance")
    # tbl_detunings = file_save.create_table(grp_detunings,"detunings",Detunings)
    # detunigns_save = tbl_detunings.row
    # for i in range(len(detunings_red)):
    #     detunigns_save["detun"] = detunings_red[i]
    #     detunigns_save.append()
    # tbl_detunings.flush()



    # num_processes = 1
    # processes = []

    # for x in range(num_processes):
    #     p = mpr.Process(target=simulation_func_red, args=(initialconds_red, 10*x,10*(x+1)))
    #     p.start()
    #     processes.append(p)

    # for p in processes: 
    #     p.join()
        
