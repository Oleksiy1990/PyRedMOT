import tables
import pandas as pd
import matplotlib.pyplot as plt 
import numpy as np 
import sys
#import matplotlib.animation as animation
import itertools

filename = "C:/Users/Oleksiy/Desktop/Code/PyRedMOT/resultsBlueMOT/starting30ms/pX17pY17pZ4grad55.hdf5"
#filename="/Users/oleksiy/Desktop/PythonCode/Transversecooling/resultsTC/pow30mWspeed550det1.hdf5"
datafile = tables.open_file(filename,"r")


radiiMM = [5,8,11,14,17]
nodes = ["radX%.iradY%.iradZ%.i/Data"%(radiusMM,radiusMM,radiusMM) for radiusMM in radiiMM]




loaded_nodes = [pd.read_hdf(filename,key=node,mode="r") for node in nodes] # we have to give it a file 
#name, the ndoe to read, and tell it to do it read-only



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

