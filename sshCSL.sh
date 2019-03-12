# setup ssh-keygen
# ssh-keygen -t rsa # hit enter until it finishes, now you will have a keygen at /Users/HHH/.ssh/id_rsa.pub
# ssh John@192.17.102.31 mkdir -p .ssh # create a .ssh folder on remote machine
# cat .ssh/id_rsa.pub | ssh John@192.17.102.31 'cat >> .ssh/authorized_keys' # append my public key to remote machine's authorized keys

# ssh to CSL desktop
# connect to VPN
ssh John@192.17.102.31
ssh -t John@192.17.102.31 "cd /home/John/riprap_project/haohang && bash" # directly ssh to a folder
# ssh -X John@192.17.102.31 # -X does not work on my mac due to missing XQuartz

# mount/unmount ssh file system
sshfs John@192.17.102.31:/home/John/riprap_project/haohang /Users/HHH/Desktop/remote
umount /Users/HHH/Desktop/remote

# transfer files to (scp [src] [dest] - 'secure copy')
scp Desktop/file.c John@192.17.102.31:/home/John/riprap_project/haohang # for files
scp -r Desktop/folder John@192.17.102.31:/home/John/riprap_project/haohang # -r for folder

# transfer files from 
scp John@192.17.102.31:/home/John/riprap_project/haohang Desktop/file.c

# Basic info
cat /etc/os-release # OS
cat /proc/cpuinfo # CPU
nvidia-smi # GPU System Management Interface
watch -n1 nvidia-smi # monitor gpu every 1 second

# Usage info
htop # or top

# Limitation
# no root priviledge
# pip --user

# ----------------------- Usage ------------------------------ #
pip list
pip show tensorflow

conda env list
source activate tf_venv
source deactivate

cat /usr/local/cuda/version.txt # CUDA version 8.0
cat /usr/local/cuda/include/cudnn.h | grep CUDNN_MAJOR -A 2 # cuDNN version 6
# Check here for compatability: https://www.tensorflow.org/install/source#tested_build_configurations
# On CSL Desktop, since the CUDA version is too old, we need downgrade:
# pip install tensorflow-gpu==1.4.0
# pip install keras=2.0.8
# other modification includes:
# mrcnn/model.py:2173 & 2197 keepdims-->keep_dims Ref: https://github.com/matterport/Mask_RCNN/issues/572


# Mask R-CNN
# Train (first cd to balloon.py, otherwise it can't find mrcnn correctly)
source activate tf_venv
export LD_LIBRARY_PATH=/usr/local/cuda-8.0/lib64
cd samples/balloon
python balloon.py train --dataset=../../datasets/balloon --weights=coco

# DeepFill
python flist_generator.py
source activate tf_venv
export LD_LIBRARY_PATH=/usr/local/cuda-8.0/lib64
python train.py


