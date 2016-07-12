import numpy as np 
import scipy as sp 
import matplotlib.pyplot as plt
import profile
import h5py
import tables #this is PyTables to use HDF5 for Python 
import sys
import os

import multiprocessing as mpr 



from scipy.integrate import odeint 
from numpy.random import random_sample, uniform
from mpl_toolkits.mplot3d import Axes3D

from pyConstants import *
from SrConstants import * 


lck = mpr.Lock()

def MOT_force(k_laser,linewidth,sat_param,velocity,B_gradient,position,detunings_list):
    delta_min = [detuning - k_laser*velocity - (muBohr*B_gradient*position/hbar) for detuning in detunings_list]
    delta_plus = [detuning + k_laser*velocity + (muBohr*B_gradient*position/hbar) for detuning in detunings_list]
    force_list = [(hbar*k_laser*linewidth*sat_param/2)*(1/(1+sat_param*(2*delta_min[q]/linewidth)**2) - \
        1/(1+sat_param*(2*delta_plus[q]/linewidth)**2)) for q in range(len(delta_min))]
    return np.sum(force_list)

def diffeqs_red(variables,t,params):
    (x,vx,y,vy,z,vz) = variables
    (laserPowerX,laserPowerY,laserPowerZ,beamWaistRadX,beamWaistRadY,beamWaistRadZ,B_gradient,detunings_list) = params


    gaussian_Xbeam = np.exp((-2*x**2)/beamWaistRadX**2)*(2*laserPowerX)/(np.pi*(beamWaistRadX**2))
    gaussian_Ybeam = np.exp((-2*y**2)/beamWaistRadY**2)*(2*laserPowerY)/(np.pi*(beamWaistRadY**2))
    gaussian_Zbeam = np.exp((-2*z**2)/beamWaistRadZ**2)*(2*laserPowerZ)/(np.pi*(beamWaistRadZ**2))


    derivs = [vx,(1/mSr88)*MOT_force(kVecRed,redGamma,gaussian_Xbeam/redIsat,vx,B_gradient,x,detunings_list),\
            vy,(1/mSr88)*MOT_force(kVecRed,redGamma,gaussian_Ybeam/redIsat,vy,B_gradient,y,detunings_list),\
            vz,(1/mSr88)*MOT_force(kVecRed,redGamma,gaussian_Zbeam/redIsat,vz,2*B_gradient,z,detunings_list) - 9.81]
            # NOTE I removed the minus sign in the z-comp of B-field 


    # Quadrupole field from: 
    # https://www2.physics.ox.ac.uk/sites/default/files/2013-01-19/minsung_pdf_16672.pdf eq. 2.40
    return derivs


def diffeqs_blue(variables,t,params):
    (x,vx,y,vy,z,vz) = variables
    (laserPowerX,laserPowerY,laserPowerZ,beamWaistRadX,beamWaistRadY,beamWaistRadZ,B_gradient,detunings_list) = params


    gaussian_Xbeam = np.exp((-2*x**2)/beamWaistRadX**2)*(2*laserPowerX)/(np.pi*(beamWaistRadX**2))
    gaussian_Ybeam = np.exp((-2*y**2)/beamWaistRadY**2)*(2*laserPowerY)/(np.pi*(beamWaistRadY**2))
    gaussian_Zbeam = np.exp((-2*z**2)/beamWaistRadZ**2)*(2*laserPowerZ)/(np.pi*(beamWaistRadZ**2))


    derivs = [vx,(1/mSr88)*MOT_force(kVecBlue,blueGamma,gaussian_Xbeam/blueIsat,vx,B_gradient,x,detunings_list),\
            vy,(1/mSr88)*MOT_force(kVecBlue,blueGamma,gaussian_Ybeam/blueIsat,vz,B_gradient,y,detunings_list),\
            vz,(1/mSr88)*MOT_force(kVecBlue,blueGamma,gaussian_Zbeam/blueIsat,vz,-2*B_gradient,z,detunings_list)]
    return derivs





def simulation_func_red(inits,start_init,end_init):
    print("Getting ready to run the solution loop")


    print("Process id: ",os.getpid())
    counter = start_init


    for init in inits[start_init:end_init]:
    

        #tbl_inits = file_save.create_table(grp_simulation,"init"+str(counter),Coords_and_vel,"Initial conditions")
        tbl_results = file_save.create_table(grp_simulation,str(counter),Coords_and_vel,"Results") 
        tbl_results.attrs.initialconds = "x: %.f, vx: %.f, y: %.f, vy: %.f, z: %.f, vz: %.f"%(init[0],init[1],init[2],init[3],init[4],init[5])
        coords_and_vel = tbl_results.row


        
        print("Solving for initial conditions %.i out of %.i"%(counter,len(inits[start_init:end_init])))
        solution_red = odeint(diffeqs_red, init, t, args=(parameters_red,),mxstep=10**8)
        print("Saving data for initial conditions %.i out of %.i"%(counter,len(inits[start_init:end_init])))

        for i in range(len(solution_red)):
            coords_and_vel["time"] = t[i]
            coords_and_vel["x_pos"] = solution_red[i,0]   
            coords_and_vel["vx"] = solution_red[i,1] 
            coords_and_vel["y_pos"] = solution_red[i,2] 
            coords_and_vel["vy"] = solution_red[i,3] 
            coords_and_vel["z_pos"] = solution_red[i,4] 
            coords_and_vel["vz"] = solution_red[i,5] 
            coords_and_vel.append()

        #we work with locks so that the file doesn't get corrupted by simultaneous saving
        lck.acquire()
        tbl_results.flush()
        lck.release()
        


        del tbl_results
        del solution_red
        print("Initial conditions %.i out of %.i are done"%(counter,len(initialconds_red)))
        counter += 1

if __name__ == '__main__': 


    num_redCombLines = 50
    detunings_red = np.linspace(-2*np.pi*200*10**3,-2*np.pi*5*10**6,num_redCombLines)


    # This is the simulation for red, let's disregard the blue for now

    # bluePowerX = bluePowerY = bluePowerZ = 20*10**-3
    # blueRadX = blueRadY = blueRadZ = 10*10**-3
    # blueGradient = 0.50 # T/m 

    redPowerXtotal = redPowerYtotal = redPowerZtotal = 10*10**-3
    redPowerX = redPowerY = redPowerZ = redPowerXtotal/num_redCombLines
    redRadX = redRadY = redRadZ = 3*10**-3
    redGradient = 1.15/100 # T/m 

    #parameters_blue = [bluePowerX,bluePowerY,bluePowerZ,blueRadX,blueRadY,blueRadZ,blueGradient,detunings_blue]
    parameters_red = [redPowerX,redPowerY,redPowerZ,redRadX,redRadY,redRadZ,redGradient,detunings_red]


    tStop = 0.5
    t = np.linspace(0., tStop, 10**5)

    position_pts_xy = 5 # number of random starting positions to check
    position_pts_z = 7
    position_rad = redRadX*np.sqrt(3/2) # distance from the center in which we allow the atom to start
                                        # we assume that if the atom is outside of e**-3 of the beam
                                        # it cannot be captured


    ## Let's leave the random positions for now 
    # pos_cube_x = uniform(low=0,high=position_rad,size=position_pts)
    # pos_cube_y = uniform(low=0,high=position_rad,size=position_pts)
    # pos_cube_z = uniform(low=-position_rad,high=position_rad,size=position_pts)
    #pos_cube = np.stack((pos_cube_x,pos_cube_y,pos_cube_z),axis=-1)



    # Making a grid or starting position coordinates
    pos_cube_x = np.linspace(0,position_rad,position_pts_xy)
    pos_cube_y = np.linspace(0,position_rad,position_pts_xy)
    pos_cube_z = np.linspace(-position_rad,position_rad,position_pts_z)


    xx,yy,zz = np.meshgrid(pos_cube_x,pos_cube_y,pos_cube_z,indexing="ij")
    xxf = xx.flatten()
    yyf = yy.flatten()
    zzf = zz.flatten()
    total_positions = np.column_stack((xxf,yyf,zzf))


    # fig = plt.figure()
    # ax = fig.add_subplot(111, projection='3d')
    # scat = ax.scatter(xx,yy,zz)
    # ax.set_xlabel('X Label')
    # ax.set_ylabel('Y Label')
    # ax.set_zlabel('Z Label')

    # plt.show()

    # sys.exit(0)

    vel_pts = 5
    #vel_min = 0.5
    vel_rad = 3 # m/s
    vel_cube_x = np.linspace(-vel_rad,vel_rad,vel_pts)
    vel_cube_y = np.linspace(-vel_rad,vel_rad,vel_pts)
    vel_cube_z = np.linspace(-vel_rad,vel_rad,vel_pts)


    vx,vy,vz = np.meshgrid(vel_cube_x,vel_cube_y,vel_cube_z,indexing="ij")
    vxf = vx.flatten()
    vyf = vy.flatten()
    vzf = vz.flatten()
    total_velocities = np.column_stack((vxf,vyf,vzf))
    total_velocities_ball = np.array([q for q in total_velocities if (q[0]**2+q[1]**2+q[2]**2) <= vel_rad**2])


    # print(total_velocities_ball)
    # fig = plt.figure()
    # ax = fig.add_subplot(111, projection='3d')
    # #scat = ax.scatter(vx,vy,vz)
    # scat = ax.scatter(total_velocities_ball[:,0],total_velocities_ball[:,1],total_velocities_ball[:,2])
    # ax.set_xlabel('X Label')
    # ax.set_ylabel('Y Label')
    # ax.set_zlabel('Z Label')

    # plt.show()

    # sys.exit(0)




    ## Let's leave the random sampling of velocities for now
    # vel_pts = 50 # number of random starting velocities to check
    # vel_rad = 5 # the max speed which will be accepted
    # vel_cube = uniform(low=-vel_rad,high=vel_rad,size=(vel_pts,1,3))
    # vel_ball = np.array([q for q in vel_cube if (q[:,0]**2+q[:,1]**2+q[:,2]**2) <= vel_rad**2])

    # fig = plt.figure()
    # ax = fig.add_subplot(111, projection='3d')
    # ax.scatter(pos_ball[:,:,0], pos_ball[:,:,1], pos_ball[:,:,2])
    # ax.set_xlabel('X Label')
    # ax.set_ylabel('Y Label')
    # ax.set_zlabel('Z Label')

    # plt.show()

    #initial_red = [(q[:,0],q[:,1],q[:,2]) for q in pos_ball]
    #initial_red = np.array([(q[:,0],v[:,0],q[:,1],v[:,1],q[:,2],v[:,2]) for q in pos_ball for v in vel_ball])


    initialconds_red = [(q[0],v[0],q[1],v[0],q[2],v[0]) for q in total_positions for v in total_velocities_ball]

    # print(len(initialconds_red))
    # sys.exit(0)

    class Coords_and_vel(tables.IsDescription):
        time = tables.Float64Col()
        x_pos = tables.Float64Col()   
        vx = tables.Float64Col()
        y_pos = tables.Float64Col()
        vy = tables.Float64Col()
        z_pos = tables.Float64Col()
        vz = tables.Float64Col()
        # idnumber  = Int64Col()      # Signed 64-bit integer
        # ADCcount  = UInt16Col()     # Unsigned short integer
        # TDCcount  = UInt8Col()      # unsigned byte
        # grid_i    = Int32Col()      # 32-bit integer
        # grid_j    = Int32Col()      # 32-bit integer
        # pressure  = Float32Col()    # float  (single-precision)
        # energy    = Float64Col()    # double (double-precision)

    class Detunings(tables.IsDescription):
        detun = tables.Float64Col()
        

    file_save = tables.open_file("C:/Users/Oleksiy/Desktop/SimulationResults/redMOT/ProcessTest_UnifPowerUnifWaistComb%.i.hdf5"%num_redCombLines,mode="a",title= "Red MOT simulation, %.i comb"%num_redCombLines)
    grp_grad = file_save.create_group("/","grad%.3f"%(redGradient*100),title="Red gradient: %.3f G/cm"%(redGradient*100))

    grp_simulation = file_save.create_group(grp_grad,str(int(redPowerXtotal*10**3))+"_"+str(int(redRadX*10**3)),\
        title="Power: " + str(int(redPowerXtotal*10**3))+ " mW, beam rad.:" +str(int(redRadX*10**3))+" mm")

    grp_detunings = file_save.create_group(grp_simulation,"detuningsHz",\
        title="Values of freqs of the comb lines, in Hz, red-detuned from resonance")
    tbl_detunings = file_save.create_table(grp_detunings,"detunings",Detunings)
    detunigns_save = tbl_detunings.row
    for i in range(len(detunings_red)):
        detunigns_save["detun"] = detunings_red[i]
        detunigns_save.append()
    tbl_detunings.flush()



    num_processes = 1
    processes = []

    for x in range(num_processes):
        p = mpr.Process(target=simulation_func_red, args=(initialconds_red, 10*x,10*(x+1)))
        p.start()
        processes.append(p)

    for p in processes: 
        p.join()
        

    # file_save.close()

#---------------------------------------------------------

