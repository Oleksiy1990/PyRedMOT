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

def MOT_force(k_laser,linewidth,sat_param,velocity,B_gradient,position,detunings_list):
    delta_min = [detuning - k_laser*velocity - (muBohr*B_gradient*position/hbar) for detuning in detunings_list]
    delta_plus = [detuning + k_laser*velocity + (muBohr*B_gradient*position/hbar) for detuning in detunings_list]
    force_list = [(hbar*k_laser*linewidth*sat_param/2)*(1/(1+sat_param*(2*delta_min[q]/linewidth)**2) - \
        1/(1+sat_param*(2*delta_plus[q]/linewidth)**2)) for q in range(len(delta_min))]
    return np.sum(force_list)



def diffeqs_blue(variables,t,params): #We disregard gravity
    (x,vx,y,vy,z,vz) = variables
    (laserPowerX,laserPowerY,laserPowerZ,beamWaistRadX,beamWaistRadY,beamWaistRadZ,B_gradient,detunings_list) = params


    gaussian_Xbeam = np.exp((-2*x**2)/beamWaistRadX**2)*(2*laserPowerX)/(np.pi*(beamWaistRadX**2))
    gaussian_Ybeam = np.exp((-2*y**2)/beamWaistRadY**2)*(2*laserPowerY)/(np.pi*(beamWaistRadY**2))
    gaussian_Zbeam = np.exp((-2*z**2)/beamWaistRadZ**2)*(2*laserPowerZ)/(np.pi*(beamWaistRadZ**2))


    derivs = [vx,(1/mSr88)*MOT_force(kVecBlue,blueGamma,gaussian_Xbeam/blueIsat,vx,B_gradient,x,detunings_list),\
            vy,(1/mSr88)*MOT_force(kVecBlue,blueGamma,gaussian_Ybeam/blueIsat,vy,B_gradient,y,detunings_list),\
            vz,(1/mSr88)*MOT_force(kVecBlue,blueGamma,gaussian_Zbeam/blueIsat,vz,-2*B_gradient,z,detunings_list)]

    # Quadrupole field from: 
    # https://www2.physics.ox.ac.uk/sites/default/files/2013-01-19/minsung_pdf_16672.pdf eq. 2.40
    return derivs



# Simulation parameters

detunings_blue = [-blueGamma]
#detunings_red = [-2*np.pi*200*10**3]

# num_redCombLines = 50
# detunings_red = np.linspace(-2*np.pi*200*10**3,-2*np.pi*5*10**6,num_redCombLines)

bluePowerX = bluePowerY = 10*10**-3 
bluePowerZ = 10*10**-3
blueRadX = blueRadY = 10*10**-3
blueRadZ = 10*10**-3
blueGradient = 0.55 # T/m = 55 G/cm



parameters_blue = [bluePowerX,bluePowerY,bluePowerZ,blueRadX,blueRadY,blueRadZ,blueGradient,detunings_blue]



tStop = 0.5
t = np.linspace(0., tStop, 10**5)

initz = np.linspace(-12e-3,12e-3,7)
init_rad = 0.05
init_angle = 30*(np.pi/180) #not to change
x_shifts = np.linspace(0,13.86e-3,5)
y_shifts = np.linspace(0,24e-3,5)
initx = -init_rad*np.sin(init_angle)+x_shifts
inity = -init_rad*np.cos(init_angle)+y_shifts

#print(initx,"\n",inity,"\n",initz)

init_speed_xy = 30 #m/s
vel_angle = np.linspace(25,35,5) #deg
initvx = init_speed_xy*np.sin(vel_angle*np.pi/180)
initvy = init_speed_xy*np.cos(vel_angle*np.pi/180)
initvz = np.linspace(0,2,4) 

inits_pos_1 = np.array(list(itertools.product(initx,[inity[0]],initz)))
inits_pos_2 = np.array(list(itertools.product([initx[0]],inity[1:],initz)))
inits_pos = np.concatenate((inits_pos_1,inits_pos_2))
print(inits_pos.shape)
inits_vel = np.array([[init_speed_xy*np.sin(angle*np.pi/180),init_speed_xy*np.cos(angle*np.pi/180),vz] for angle in vel_angle for vz in initvz])
print(inits_vel.shape)

inits = np.array([[p[0],v[0],p[1],v[1],p[2],v[2]] for p in inits_pos for v in inits_vel])
print(inits.shape)



final_pos_captured = []
final_pos_lost = []
weird = []

#inits = (initx,initvx,inity,initvy,initz,initvz)
for num,ic in enumerate(inits):
	print("Solving ",num)
	solution_blue = odeint(diffeqs_blue, ic, t, args=(parameters_blue,),mxstep=10**8)
	fp = np.sqrt(solution_blue[-1,0]**2+solution_blue[-1,2]**2+solution_blue[-1,4]**2)
	if fp<=1e-4:
		final_pos_captured.append(fp)
	elif fp>1e-4:
		final_pos_lost.append(fp)
	else:
		weird.append(fp)

print(len(final_pos_captured))
print(len(final_pos_lost))
print(len(weird))

sys.exit(0)

position_res = np.sqrt(solution_blue[:,0]**2+solution_blue[:,2]**2+solution_blue[:,4]**2)

# Plot results
fig = plt.figure(1, figsize=(8,8))

# Plot theta as a function of time
ax1 = fig.add_subplot(111)
ax1.scatter(solution_blue[:,0],solution_blue[:,2])
ax1.set_xlabel('x')
ax1.set_ylabel('y')
plt.show()


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

