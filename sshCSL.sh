# setup ssh-keygen
# ssh-keygen -t rsa # hit enter until it finishes, now you will have a keygen at /Users/HHH/.ssh/id_rsa.pub
# ssh John@192.17.102.31 mkdir -p .ssh # create a .ssh folder on remote machine
# cat .ssh/id_rsa.pub | ssh John@192.17.102.31 'cat >> .ssh/authorized_keys' # append my public key to remote machine's authorized keys

# ssh to CSL desktop
# connect to VPN
ssh John@192.17.102.31
ssh -t John@192.17.102.31 "cd /home/John/riprap_project/haohang && bash" # directly ssh to a folder
# ssh -Y John@192.17.102.31 # -Y with GUI, first start XQuartz on mac before type in terminal

# check running process
tmux attach
# detach
press ctrl+B  and press D

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
# pip install tensorflow==1.5.0
# pip install tensorflow-gpu==1.4.0
# pip install keras=2.0.8
# other modification includes:
# mrcnn/model.py:2173 & 2197 keepdims-->keep_dims Ref: https://github.com/matterport/Mask_RCNN/issues/572

# Teamviewer: https://www.techinfected.net/2016/10/install-use-teamviewer-on-ubuntu.html
wget https://download.teamviewer.com/download/teamviewer_i386.deb
sudo dpkg -i teamviewer*.deb
sudo apt-get -f install

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

# VisualSfM: http://ccwu.me/vsfm/
# http://www.cs.cmu.edu/~reconstruction/basic_workflow.html

# Installation Guide
# Download SiftGPU, pba, PMVS, CMVS, Graclus, compile them and copy the compiled library to /vsfm/bin/
# Ref: https://blog.csdn.net/binghunwangcong/article/details/79866093
# Ref: https://gist.github.com/lvisintini/e07abae48f099b913f9cf1c1f0fe43ba
# Ref: http://www.10flow.com/2012/08/15/building-visualsfm-on-ubuntu-12-04-precise-pangolin-desktop-64-bit/
# Ref: https://www.cnblogs.com/gaoxiang12/p/5149067.html

# Dependencies
sudo apt install \
	libglib2.0-dev \
	libglib2.0-0=2.48.2-0ubuntu1 \
    libgtk2.0-dev \
    libglew-dev \
    libdevil-dev \
    libboost-all-dev \
    libatlas-cpp-0.6-dev \
    libatlas-dev \
    imagemagick \
    libatlas3-base \
    libcminpack-dev \
    libgfortran3 \
    libmetis-edf-dev \
    libparmetis-dev \
    freeglut3-dev \
    libgsl-dev \
    glew-utils \
    liblapack-dev

# Download VisualSFM: http://ccwu.me/vsfm/
wget http://ccwu.me/vsfm/download/VisualSFM_linux_64bit.zip
unzip VisualSFM_linux_64bit.zip
rm VisualSFM_linux_64bit.zip

# Download SiftGPU: https://github.com/pitzer/SiftGPU
cd to /vsfm
wget https://github.com/pitzer/SiftGPU/archive/master.zip
unzip master.zip
rm master.zip
mv SiftGPU-master SiftGPU

cd to /vsfm/SiftGPU
make
cp /vsfm/SiftGPU/bin/libsiftgpu.so to /vsfm/bin # put the compiled library libsiftgpu.so to /vsfm/bin

# Download Multicore Bundle Adjustment ('pba'): http://grail.cs.washington.edu/projects/mcba/
cd to /vsfm
wget http://grail.cs.washington.edu/projects/mcba/pba_v1.0.5.zip
unzip pba_v1.0.5.zip
rm pba_v1.0.5.zip

cd to /vsfm/pba
# Try using the following commands if the make command failed
# echo -e "#include <stdlib.h>\n$(cat ~/vsfm/pba/src/pba/SparseBundleCU.h)" > ~/vsfm/pba/src/pba/SparseBundleCU.h
# echo -e "#include <stdlib.h>\n$(cat ~/vsfm/pba/src/pba/pba.h)" > ~/vsfm/pba/src/pba/pba.h
make
cp /vsfm/pba/bin/libpba.so to /vsfm/bin/

# Download PMVS 
# http://www.di.ens.fr/pmvs/documentation.html
cd to vsfm
wget http://www.di.ens.fr/pmvs/pmvs-2.tar.gz
tar xvzf pmvs-2.tar.gz
rm pmvs-2.tar.gz

cd to /vsfm/pmvs-2/program/main/
cp ~/vsfm/pmvs-2/program/main/mylapack.o ~/vsfm/pmvs-2/program/main/mylapack.o.backup
make clean
cp ~/vsfm/pmvs-2/program/main/mylapack.o.backup ~/vsfm/pmvs-2/program/main/mylapack.o
make depend
make

# problem missing -llapack (solved)
# apt install liblapack-dev
# https://blog.csdn.net/binghunwangcong/article/details/79866093

# Download Graclus1.2: http://www.cs.utexmcas.edu/users/dml/Software/graclus.html
cd to /vsfm
wget http://www.cs.utexas.edu/users/dml/Software/graclus1.2.tar.gz
tar xvzf graclus1.2.tar.gz
rm graclus1.2.tar.gz
sed -i 's/COPTIONS = -DNUMBITS=32/COPTIONS = -DNUMBITS=64/' ~/vsfm/graclus1.2/Makefile.in
cd to /vsfm/graclus1.2/
make

# Download CMVS: http://www.di.ens.fr/cmvs/documentation.html
cd to /vsfm
wget http://www.di.ens.fr/cmvs/cmvs-fix2.tar.gz
tar xvzf cmvs-fix2.tar.gz
rm cmvs-fix2.tar.gz

cp /vsfm/pmvs-2/program/main/mylapack.o /vsfm/cmvs/program/main/

echo -e "#include <vector>\n#include <numeric>\n$(cat ~/vsfm/cmvs/program/base/cmvs/bundle.cc)" > ./vsfm/cmvs/program/base/cmvs/bundle.cc # add two lines in bundle.cc
echo -e "#include <stdlib.h>\n$(cat ~/vsfm/cmvs/program/main/genOption.cc)" > ./vsfm/cmvs/program/main/genOption.cc # add one line in genOption.cc
sed -e '/Your INCLUDE path*/ s/^#*/#/' -i ./vsfm/cmvs/program/main/Makefile
sed -e '/Your metis directory*/ s/^#*/#/' -i ./vsfm/cmvs/program/main/Makefile
sed -e '/Your LDLIBRARY path*/ s/^#*/#/' -i ./vsfm/cmvs/program/main/Makefile

sed -i "s:YOUR_INCLUDE_METIS_PATH =*:YOUR_INCLUDE_METIS_PATH = -I$HOME/vsfm/graclus1.2/metisLib:" ./vsfm/cmvs/program/main/Makefile # should ensure you have the right path to graclus1.2
sed -i "s:YOUR_LDLIB_PATH =*:YOUR_LDLIB_PATH = -L$HOME/vsfm/graclus1.2:" ./vsfm/cmvs/program/main/Makefile

cd to /vsfm/cmvs/program/main
make

cp ~/vsfm/cmvs/program/main/cmvs ~/vsfm/bin
cp ~/vsfm/cmvs/program/main/pmvs2 ~/vsfm/bin
cp ~/vsfm/cmvs/program/main/genOption ~/vsfm/bin

cd to /vsfm # now /vsfm/bin/ should have cmvs  genOption  libpba.so  libsiftgpu.so  pmvs2
make # now you got a VisualSfM executable

# missing libgtk2.0-dev
# https://askubuntu.com/questions/884413/unable-to-install-libgtk2-0-dev-on-ubuntu-16-04
sudo apt-get install libglib2.0-dev libglib2.0-0=2.48.2-0ubuntu1
sudo apt-get install libgtk2.0-dev
# missing liblapack-dev
sudo apt install liblapack-dev

# Teamviewer Install
wget https://download.teamviewer.com/download/teamviewer_i386.deb
sudo dpkg -i teamviewer*.deb
sudo apt-get -f install

# Teamviewer Use
systemctl status teamviewerd # check if teamviewer is running
sudo systemctl start teamviewerd # manually start the daemon
sudo systemctl enable teamviewerd # autostart at boot time

teamviewer --passwd 3212 # set password at 1st run, password is password
teamviewer -info # show the ID

