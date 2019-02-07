# ================================ Argv ========================================
import argparse
parser = argparse.ArgumentParser()
parser.add_argument("-f", "--filename", help="file name")
args = parser.parse_args()
filename = args.filename
# ============================== Compressor ====================================
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# Read from .obj file
v_coords = []
with open(filename, 'r') as f:
    lines = f.read().splitlines()[2:]
    for l in range(len(lines)):
        tokens = lines[l].split(' ')
        if (tokens[0] == 's'):
            break
        for t in tokens[1:]:
            v_coords.append(float(t))

# Compress to [0,1] range
point_cloud = np.asarray(v_coords, dtype = float).reshape((-1,3))
coords_min = np.amin(point_cloud, axis = 0)
coords_max = np.amax(point_cloud, axis = 0)
point_cloud_compressed = (point_cloud - coords_min) / (coords_max - coords_min) # broadcasting
# point_cloud_compressed is the N x 3 matrix of point cloud coordinates

# Rasterization
x_coords = point_cloud_compressed[:,0]
y_coords = point_cloud_compressed[:,1]
z_coords = point_cloud_compressed[:,2]

# Plot
PLOT = False
if PLOT:
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.scatter(x_coords, y_coords, z_coords, zdir = 'z', c = 'r')
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    # plt.show()
    # plt.savefig("test.png")
