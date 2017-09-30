import tables
import pandas as pd
import matplotlib.pyplot as plt 
import numpy as np 
import sys
#import matplotlib.animation as animation
import itertools

#filename = "C:/Users/Oleksiy/Desktop/Code/PyRedMOT/resultsBlueMOT/starting30ms/pX17pY17pZ4grad55.hdf5"

directory = "/Users/oleksiy/Desktop/PythonCode/PyRedMOT/resultsBlueMOT/"

simNumber = 8
filepX7pY7pZ7 ="sim%.ipX7pY7pZ7grad24.hdf5"%simNumber
filepX8pY8pZ5 ="sim%.ipX8pY8pZ5grad24.hdf5"%simNumber
filepX10pY10pZ2 ="sim%.ipX10pY10pZ2grad24.hdf5"%simNumber

beamRadiiMM = [1,3,5,7,9,11,13,15,17]
# each entry in the following lists corresponds to different beam radii, as given in the variable beamRadiiMM
dfsX7Y7Z7 = [pd.read_hdf(directory+filepX7pY7pZ7,"radX%.iradY%.iradZ%.i/Data"%(rad,rad,rad)) for rad in beamRadiiMM]
dfsX8Y8Z5 = [pd.read_hdf(directory+filepX8pY8pZ5,"radX%.iradY%.iradZ%.i/Data"%(rad,rad,rad)) for rad in beamRadiiMM]
dfsX10Y10Z2 = [pd.read_hdf(directory+filepX10pY10pZ2,"radX%.iradY%.iradZ%.i/Data"%(rad,rad,rad)) for rad in beamRadiiMM]  

#distances are in m and velocities are in m/s, so all SI units
stoppedX7Y7Z7 = [q[(q.x_final)**2 + (q.y_final)**2 + (q.z_final)**2 < 0.001**2] for q in dfsX7Y7Z7]
stoppedX8Y8Z5 = [q[(q.x_final)**2 + (q.y_final)**2 + (q.z_final)**2 < 0.001**2] for q in dfsX8Y8Z5]
stoppedX10Y10Z2 = [q[(q.x_final)**2 + (q.y_final)**2 + (q.z_final)**2 < 0.001**2] for q in dfsX10Y10Z2]

num_stoppedX7Y7Z7 = [len(q) for q in stoppedX7Y7Z7]
num_stoppedX8Y8Z5 = [len(q) for q in stoppedX8Y8Z5]
num_stoppedX10Y10Z2 = [len(q) for q in stoppedX10Y10Z2]

plt.scatter(beamRadiiMM,num_stoppedX7Y7Z7,label="Power X7 Y7 Z7 mW",c="g",marker="o")
plt.scatter(beamRadiiMM,num_stoppedX8Y8Z5,label="Power X8 Y8 Z5 mW",c="b",marker="v")
plt.scatter(beamRadiiMM,num_stoppedX10Y10Z2,label="Power X10 Y10 Z2 mW",c="r",marker="s")
plt.xlabel("Beam radius")
plt.ylabel("Number of initial conditions captured atoms")
plt.legend(loc='upper center',bbox_to_anchor=(0.5, 1.14))
plt.show()

print(num_stoppedX7Y7Z7)
print(num_stoppedX8Y8Z5)
print(num_stoppedX10Y10Z2)

sys.exit(0)



"""
We load the HDF5 which has beed read before as a DataFrame (using Pandas)

We then check condition by condition to narrow down and select the data that we like 
"""
def distance_func(vec):
    return np.sqrt(vec[0]**2+vec[1]**2+vec[2]**2)

final_positions = []
captured_number = []

for index,u in enumerate(loaded_nodes):

    final_positions.append(pd.DataFrame(u,columns=["x_final","y_final","z_final"]))

for q in final_positions:
    captured = q[abs(q["x_final"])<1e-4]
    captured = captured[abs(captured["y_final"])<1e-4]
    captured = captured[abs(captured["z_final"])<1e-4]
    print(captured)
    captured_number.append(len(captured))
    
print(captured_number)

sys.exit(0)

final_positions = np.array(final_positions)


sys.exit(0)
    

    # success_data = full_data
    # #success_data = success_data[success_data['init_x'] < 2.1e-3]
    # #success_data = success_data[success_data['init_y'] < 1e-3]
    # #success_data = success_data[success_data['init_vy'] < 1e-3]
    # #success_data = success_data[success_data['init_vx'] > 2]
    # #success_data = success_data[success_data['init_vx'] < 10]
    # success_data = success_data[success_data['final_speed_xy'] <= 1.5]
    # #print(np.sqrt(success_data['init_vx']**2+success_data['init_vx']**2))
    
    
    # fail_data = full_data
    # #fail_data = fail_data[fail_data['init_x'] < 2.1e-3]
    # #success_data = success_data[success_data['init_y'] < 1e-3]
    # #success_data = success_data[success_data['init_vy'] < 1e-3]
    # #success_data = success_data[success_data['init_vx'] > 2]
    # #success_data = success_data[success_data['init_vx'] < 10]
    # fail_data = fail_data[fail_data['final_speed_xy'] > 1.5]
    # #print(success_data)
    
    # fig = plt.figure()
    # ax1 = fig.add_subplot(211)
    
    
    # ax1.scatter(np.sqrt(success_data["init_vx"]**2 + success_data["init_vy"]**2),success_data["final_speed_xy"],c="g")
    # ax1.set_title("a%.i b%.i"%ab_combinations[index])
    # ax1.set_xlabel("Init speed [m/s]")
    # ax1.set_ylabel("Final speed [m/s]")
    
    # ax2 = fig.add_subplot(212)
    
    # ax2.scatter(np.sqrt(fail_data["init_vx"]**2 + fail_data["init_vy"]**2),fail_data["final_speed_xy"],c="r")
    # ax2.set_title("a%.i b%.i"%ab_combinations[index])
    # ax2.set_xlabel("Init speed [m/s]")
    # ax2.set_ylabel("Final speed [m/s]")
    
    # plt.show()
    
    # print("Captured: %.i"%len(success_data))
    # print("Failed: %.i"%len(fail_data))
    # print("Ratio: %.4f"%(len(success_data)/len(fail_data)))
    


"""


Notes: 
"index" apparently means rows and "columns" are as it sounds in Pandas

"""

