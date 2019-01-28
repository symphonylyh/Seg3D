#!/bin/bash

# Add meshlabserver directory and libraries to OS PATH
export MESHLABSERVER_PATH=/Applications/meshlab.app/Contents/MacOS # meshlabserver PATH
export DYLD_FRAMEWORK_PATH=/Applications/meshlab.app/Contents/Frameworks # meshlabserver's dynamic library path

/Applications/meshlab.app/Contents/MacOS/meshlabserver -d filters.txt
# /Applications/meshlab.app/Contents/MacOS/meshlabserver -i test.ply -o test.off -w outproj.mlp -x -s test.mlx