# Step 1: Install Blender (v2.79)
#
# Step 2: Download and Install rock generator add-on 'add_mesh_rocks.zip'
# 	Tutorial: https://blender.stackexchange.com/questions/118887/where-is-the-blender-rock-generator-addon
# 	File-User Preferences-Add-ons-Install Add-on from file
# 	The function will then be enabled under Add-Mesh-Rock Generator
# 	Doc: https://archive.blender.org/wiki/index.php/Extensions:2.6/Py/Scripts/Add_Mesh/Rock_Generator/
#
# Step 3: Run this script to Generate and Save rock models
# https://www.youtube.com/watch?v=K0yb4sZ7B4g
# https://blender.stackexchange.com/questions/33755/batch-exporting-of-multiple-objects-into-separate-stl-files
#
# Step 4: Run MeshLab to get the point cloud data (actually unnecessary, directly parse the .obj file to get point cloud)

# ================================ Argv ========================================
# https://blender.stackexchange.com/questions/6817/how-to-pass-command-line-arguments-to-a-blender-python-script
import sys
argv = sys.argv
argv = argv[argv.index("--") + 1:]  # get all args after "--"
N = int(argv[0]) # No. of rocks
PATH = argv[1] # Destination path

# ============================== Generator =====================================
import bpy # Python API of Blender
import os # mkdir

# Define Directory (Create Directory If Necessary)
dir = bpy.path.abspath(PATH)
# if not os.path.exists(dir): os.makedirs(dir)

# Remove 'camera', 'lamp' and 'cube'
scene = bpy.context.scene
for ob in scene.objects:
    if ob.type == 'CAMERA' or ob.type == 'LAMP' or (ob.type == 'MESH' and ob.name.startswith("Cube")):
        ob.select = True
    else:
        ob.select = False
bpy.ops.object.delete()

# Generate rocks
bpy.ops.mesh.rocks(num_of_rocks = N) # in Blender GUI, floating over the toggle, the corresponding variable name will pop up as "MEST_OT.num_of_rocks"

# Export mesh (.obj)
for object in bpy.context.scene.objects:
    # deselect all meshes
    bpy.ops.object.select_all(action='DESELECT')
    # select the object
    object.select = True
    # export object with its name as file name
    filename = str((dir + object.name + '.obj'))
    #bpy.context.active_object = object
    bpy.ops.export_scene.obj(filepath = filename, use_selection = True, use_edges = False, use_normals = False, use_uvs = False, use_materials = False, use_triangles = True, use_blen_objects = False)
    # Doc: https://docs.blender.org/api/2.79/bpy.ops.export_scene.html

    # Export mesh (.stl) or (.ply) .stl has batch mode, .ply will merge all rocks into one...(damn .ply)
    # bpy.ops.export_mesh.stl(filepath = '/Users/HHH/Desktop/rocks/', ascii = True, batch_mode = 'OBJECT') # stl has batch mode
    # Doc: https://docs.blender.org/api/2.79/bpy.ops.export_mesh.html
