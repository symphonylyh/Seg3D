# Step 1 Load Imageset
# File--Open+ Multi Images

# Step 2 Detect SIFT features
# SfM--Pairwise Matching--Compute Missing Match

# Step 3 Sparse Reconstruction
# SfM--Reconstruct Sparse

# Step 4 Dense Reconstruction
# SfM--Reconstruct Dense

# Optional Steps for Scale Calibration
# Option 1:
# SfM--More Fucntions--GCP-Based Transform
# Ctrl+Click points: choose points and specify the real-world coordinates
# Shift+Click GCP label: compute transformation matrix
# X' = S * R * X + T.
# X is the user-defined coordinates system according to real-world coordinates
# X' is the current SfM coordinates system
# S is scale factor, R is rotation matrix, T is translation vector
# so current SfM scale = S * use-defined scale
# i.e. we need to apply 1/S to rescale the point cloud data
# OR
# MeshLab measurement
# Option 2:
# Use a fixed two-camera frame (like Kinect, distance between two cameras is known)
# Get the SfM coordinates of two cameras and scale with the measured distance

# name.nvm.cmvs/00/models/options-0000.ply is the geometry file
