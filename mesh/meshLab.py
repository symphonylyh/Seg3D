import os 
import meshlabxml as mlx 

# Add meshlabserver directory and libraries to OS PATH
# MESHLABSERVER_PATH = '/Applications/meshlab.app/Contents/MacOS' # meshlabserver PATH
# DYLD_FRAMEWORK_PATH = '/Applications/meshlab.app/Contents/Frameworks' # meshlabserver's dynamic library path
# os.environ['PATH'] += os.pathsep + MESHLABSERVER_PATH + os.pathsep + DYLD_FRAMEWORK_PATH
os.environ['MESHLABSERVER_PATH'] = "/Applications/meshlab.app/Contents/MacOS"
os.environ['DYLD_FRAMEWORK_PATH'] = "/Applications/meshlab.app/Contents/Frameworks"

os.system('/Applications/meshlab.app/Contents/MacOS/meshlabserver -h')

# orange_cube = mlx.FilterScript(file_out='orange_cube.ply', ml_version='2016.12')
# mlx.create.cube(orange_cube, size=[3.0, 4.0, 5.0], center=True, color='orange')
# mlx.transform.rotate(orange_cube, axis='x', angle=45)
# mlx.transform.rotate(orange_cube, axis='y', angle=45)
# mlx.transform.translate(orange_cube, value=[0, 5.0, 0])
# orange_cube.run_script()

# script = 'orange_cube.mlx' # script file 
# model = 'orange_cube.ply' # output file 
# log = 'orange_cube_log.txt' # log file 

# mlx.begin(script=script) # Start writing the script to the script file 
# mlx.create.cube(script=script, size=[3.0, 3.0, 2.0], center=True, color='orange') 
# mlx.transform.rotate(script=script, axis='x', angle=45) 
# mlx.transform.translate(script=script, value=[5.0, 0, 0]) 
# mlx.end(script=script) # Finish writing the script to the script file 

# mlx.run(script=script, log=log, file_out=model) # Run the script using meshlabserver and output the result 
# mlx.util.delete_all('TEMP3D*') # Delete temp files 