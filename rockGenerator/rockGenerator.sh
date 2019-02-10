#!/bin/bash
# Usage: sh rockGenerator.sh 100 // to generate 100 rocks

# Configuration
  # 1. have rockGenerator.sh & rockGenerator.py under the same folder
  # 2. cd to this folder
# Output
  # results will be under a new folder called 'samples' under the same folder

export BLENDER_PATH=/Applications/Blender/blender.app/Contents/MacOS

# First argv $1 = No. of rocks
[[ -d "./samples_$1" ]] || mkdir "./samples_$1"

/Applications/Blender/blender.app/Contents/MacOS/blender --background --python ./rockGenerator.py -- $1 ./samples_$1/ # pass additional arguments
