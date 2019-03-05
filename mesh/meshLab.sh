#!/bin/bash

# Add meshlabserver directory and libraries to OS PATH
export MESHLABSERVER_PATH=/Applications/meshlab.app/Contents/MacOS # meshlabserver PATH
export DYLD_FRAMEWORK_PATH=/Applications/meshlab.app/Contents/Frameworks # meshlabserver's dynamic library path

/Applications/meshlab.app/Contents/MacOS/meshlabserver -i test.ply -o test.off -w outproj.mlp -x -s test.mlx

# MeshLab Manual Procedure for 01_06_2019 dataset (using blackboard as reference object):
# Step 1: 
# Take the option-0000.ply file from VisualSfM output (this is the raw point cloud data)

# Step 2: 
# Clean the ground by removing vertices
# Poisson Reconstruction with Depth = 9

# Step 3: 
# Measure the dimension of the board and record in a .txt file (Measure directly on point cloud is nearly impossible, so measure on mesh)
# Clean the mesh to only have the rock

# Step 4:
# Save the ground removed point cloud: Export Mesh-Untick "Binary Encoding". This will overwrite the previous option-0000.ply file
# Save the MeshLab project: Save Project-Poisson Mesh layer into 01.off file, whole project into 01.mlp file, 
# .mlp & .ply files are for documenting use, .off is for Seg3D code.
# Next time, just open .mlp in MeshLab