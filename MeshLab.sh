# Step 1: Simpify point cloud
# Filters-Cleaning and Repairing-Merge Close Vertices
# Filters-Sampling-Point Cloud Simplification-Number of Samples: 100000
# OR
# Filters-Sampling-Poisson-disk sampling-Toggle Base mesh subsampling, Number of Samples: 100000

# Step 2: Surface reconstruction
# Filters-Remeshing, Simplification and Reconstruction-Screened Poisson Surface Reconstruction

# Step 3: Clean mesh
# Filters-Selection-Select faces with edge longer than...-Delete selected faces and vertices
# Filters-Selection-Small component selection-Delete selected faces and vertices

# Step 4: Scale calibration
# Measuring tool--Selection reference object--Record length and manully calculate scale factor
# Q1. How to get the real scale of scene in world units?
# Option 1: Measure the distance between two camera locations (either
# manually or connect two cameras with know distance), then this distance
# actually becomes a "virtual" calibration ruler
# Option 2: Attach a motion sensor with camera that can record the moving distance between each photo
