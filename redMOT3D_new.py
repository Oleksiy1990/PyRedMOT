import numpy as np 
import scipy as sp 
import matplotlib.pyplot as plt
import profile

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
    force_max = hbar*k_laser*linewidth/2
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

    phi = uniform(0,np.pi)
    theta = uniform(0,2*np.pi)


    derivs = [vx,(1/mSr88)*(MOT_force(k_vector,Gamma,gaussian_Xbeam/Isat,vx,B_gradient,x,detuning)+MOT_force(k_vector,Gamma,gaussian_Xbeam/Isat,vx,B_gradient,x,detuning)*np.sin(theta)*np.cos(phi)),\
            vy,(1/mSr88)*(MOT_force(k_vector,Gamma,gaussian_Ybeam/Isat,vy,B_gradient,y,detuning)+MOT_force(k_vector,Gamma,gaussian_Xbeam/Isat,vx,B_gradient,x,detuning)*np.sin(theta)*np.sin(phi)),\
            vz,(1/mSr88)*MOT_force(k_vector,Gamma,gaussian_Zbeam/Isat,vz,2*B_gradient,z,detuning)+MOT_force(k_vector,Gamma,gaussian_Xbeam/Isat,vx,B_gradient,x,detuning)*np.cos(phi) - g_Earth]

    # Quadrupole field from: 
    # https://www2.physics.ox.ac.uk/sites/default/files/2013-01-19/minsung_pdf_16672.pdf eq. 2.40
    return derivs

# Simulation parameters



# This is the simulation for red, let's disregard the blue for now


detuning_blue = -blueGamma*2
bluePowerX = bluePowerY = bluePowerZ = 10e-3
blueRadX = blueRadY = blueRadZ = 10e-3
blueGradient = -0.45 # T/m = 45 G/cm NOTE! Due to sign conventions, this must be written as negative, otherwise equations fail


parameters_blue = [kVecBlue,blueGamma,blueIsat,bluePowerX,bluePowerY,bluePowerZ, \
                    blueRadX,blueRadY,blueRadZ,blueGradient,detuning_blue]


tStop = 0.1
t = np.linspace(0., tStop, 1e5)
inits_blue_test = (0,5,0,5,0,5)
solution_blue = odeint(diffeqs_MOT, inits_blue_test, t, args=(parameters_blue,),mxstep = 10**1)
print(solution_blue)


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
plt.plot(t,solution_blue[:,1],label="vx")
plt.plot(t,solution_blue[:,3],label="vy")
plt.plot(t,solution_blue[:,5],label="vz")
plt.legend()
plt.show()
sys.exit(0)





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
    

file_save = tables.open_file("C:/Users/Oleksiy/Desktop/SimulationResults/redMOT/UnifPowerUnifWaistComb%.i.hdf5"%num_redCombLines,mode="a",title= "Red MOT simulation, %.i comb"%num_redCombLines)
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



print("Getting ready to run the solution loop")

counter = 0

for inits in initialconds_red: 

    #tbl_inits = file_save.create_table(grp_simulation,"init"+str(counter),Coords_and_vel,"Initial conditions")
    tbl_results = file_save.create_table(grp_simulation,str(counter),Coords_and_vel,"Results")
    tbl_results.attrs.initialconds = "x: %.f, vx: %.f, y: %.f, vy: %.f, z: %.f, vz: %.f"%(inits[0],inits[1],inits[2],inits[3],inits[4],inits[5])
    coords_and_vel = tbl_results.row


    print("Solving for initial conditions %.i out of %.i"%(counter,len(initialconds_red)))
    solution_red = odeint(diffeqs_red, inits, t, args=(parameters_red,),mxstep=10**8)
    print("Saving data for initial conditions %.i out of %.i"%(counter,len(initialconds_red)))

    for i in range(len(solution_red)):
        coords_and_vel["time"] = t[i]
        coords_and_vel["x_pos"] = solution_red[i,0]   
        coords_and_vel["vx"] = solution_red[i,1] 
        coords_and_vel["y_pos"] = solution_red[i,2] 
        coords_and_vel["vy"] = solution_red[i,3] 
        coords_and_vel["z_pos"] = solution_red[i,4] 
        coords_and_vel["vz"] = solution_red[i,5] 
        coords_and_vel.append()

    tbl_results.flush()
    del tbl_results
    del solution_red
    print("Initial conditions %.i out of %.i are done"%(counter,len(initialconds_red)))
    counter += 1

file_save.close()
sys.exit(0)


#---------------------------------------------------------

coords_and_vel = tbl.row

print("Getting ready to solve")
solution_red = odeint(diffeqs_red, initialconds_red[0], t, args=(parameters_red,),mxstep=10**9)

print("Solved! Getting ready to write data")

for i in range(len(solution_red)):
    coords_and_vel["x_pos"] = solution_red[i,0]   
    coords_and_vel["vx"] = solution_red[i,1] 
    coords_and_vel["y_pos"] = solution_red[i,2] 
    coords_and_vel["vy"] = solution_red[i,3] 
    coords_and_vel["z_pos"] = solution_red[i,4] 
    coords_and_vel["vz"] = solution_red[i,5] 
    coords_and_vel.append()

tbl.flush()

sys.exit(0)




def simulation_func_red(initialconds,parameters_red,t):
    solutions = []
    for inits in initialconds_red:
        solution_red = odeint(diffeqs_red, inits, t, args=(parameters_red,))
        endresult = solution_red[-1,:]
        solutions.append(endresult)
        del solution_red
        del endresult
    return solutions

# profile.run("main()")

sols = simulation_func_red(initialconds_red,parameters_red,t)












# dset_zplus_inits[...] = initialconds_red
# dset_zplus_data[...] = sols
# file_save.close()

# f = h5py.File("resultsRed/test.hdf5","a")
# grp = f.create_group("mygroup")
# dset = grp.create_dataset("mydataset",data=sols)
# grp.attrs["power_g"] = 10
# grp.attrs["width_g"] = 5
# dset.attrs["power"] = 10
# dset.attrs["width"] = 5
# dset.attrs["dets"] = detunings_red

# np.savetxt("inits1.txt",initialconds_red)
# np.savetxt("results1.txt",sols)



# # Plot results
# fig = plt.figure(1, figsize=(8,8))

# # Plot theta as a function of time
# ax1 = fig.add_subplot(211)
# ax1.plot(t, solution_red[:,4])
# ax1.set_xlabel('time')
# ax1.set_ylabel('theta')

# # Plot omega as a function of time
# ax2 = fig.add_subplot(212)
# ax2.plot(t, solution_red[:,5])
# ax2.set_xlabel('time')
# ax2.set_ylabel('omega')


# plt.tight_layout()
# plt.show()