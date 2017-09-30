import tables
import pandas as pd
import matplotlib.pyplot as plt 
from mpl_toolkits.mplot3d import Axes3D
import numpy as np 
import sys
#import matplotlib.animation as animation
import itertools


num_lines = 100
gradient = 115


filename = "C:/Users/Oleksiy/Desktop/Code/PyRedMOT/resultsRedMOT/redMOT%.ilinesGravityOn.hdf5"%num_lines
#filename="/Users/oleksiy/Desktop/PythonCode/Transversecooling/resultsTC/pow30mWspeed550det1.hdf5"

radiiMM = [2,7,12,17]
powermW = 1


nodes = ["gradient%.iGcm/Pow%.imWredRad%.immGrad%.iGcmlines%.i"%(int(gradient),powermW,radiusMM,int(gradient),num_lines) for radiusMM in radiiMM]
data = [pd.read_hdf(filename,key=node_name,mode="r") for node_name in nodes] 

# print(type(u))
# print(u.index)
# print(u.columns)

# print(type(captured))
# print(captured)


# captured = u[abs(u["x_final"])<1e-4]

fig = plt.figure()
ax = fig.add_subplot(211, projection='3d')

for pts,c in [(data[0],"b"),(data[1],"g"),(data[2],"r"),(data[3],"k")]:
    captured = pts[abs(pts["x_final"])<5e-4]
    ax.scatter(captured["x_init"],captured["vx_init"],captured["r_init"],c=c)

ax.legend(radiiMM)
ax.set_xlabel("Initial x [m]")
ax.set_ylabel("Initial vx [m/s]")
ax.set_zlabel("Initial r in beam [m]")
ax.set_title("Captured\n Power %.i mW, %.i lines"%(powermW,num_lines))

ax = fig.add_subplot(212, projection='3d')

for pts,c in [(data[0],"b"),(data[1],"g"),(data[2],"r"),(data[3],"k")]:
    captured = pts[abs(pts["x_final"])>5e-4]
    ax.scatter(captured["x_init"],captured["vx_init"],captured["r_init"],c=c)

ax.legend(radiiMM)
ax.set_xlabel("Initial x [m]")
ax.set_ylabel("Initial vx [m/s]")
ax.set_zlabel("Initial r in beam [m]")
ax.set_title("Escaped\n Power %.i mW, %.i lines"%(powermW,num_lines))


plt.show()



sys.exit(0)




