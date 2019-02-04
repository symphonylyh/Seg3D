#!/bin/bash
# Configuration
  # 1. have rockGenerator.sh & rockGenerator.py under the same folder
  # 2. cd to this folder
# Output
  # results will be under a new folder called 'samples' under the same folder

export BLENDER_PATH=/Applications/Blender/blender.app/Contents/MacOS

N=5 # No. of rocks
[[ -d "./samples_$N" ]] || mkdir "./samples_$N"

/Applications/Blender/blender.app/Contents/MacOS/blender --background --python ./rockGenerator.py -- $N ./samples_$N/ # pass additional arguments
