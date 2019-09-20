#!/bin/sh

# Run docker with current directory mounted as /data
# Usage: 
#   run_hb_lawson_docker.sh bash : run bash inside docker.
#   run_hb_lawson_docker.sh seg T1w.nii.gz my_seg: open segmentation program
#   run_hb_lawson_docker.sh to_nifti T1w.nii.gz my_seg_hb_lawson.csv: convert csv file to nifti file

if [ N`uname` = "NDarwin" ]; then
    IP=$(ifconfig en0 | grep inet | awk '$1=="inet" {print $2}') # if IP is wrong, correct this line
    echo "IP for X-forwarding: $IP"
    DISP=${IP}:0
else
    DISP=$DISPLAY
fi


USAGE="Usage:\n
    $0 bash: run bash inside docker\n
    $0 seg T1w.nii.gz my_seg: open segmentation program\n
    $0 to_nifti T1w.nii.gz my_seg_hb_lawson.csv: convert csv file to nifti file\n
"

if [ $# -eq 0 ]; then
    echo $USAGE
    exit 1
fi

if [ $1 = "seg" ]; then
    CMD="hb_segmentation_lawson.py"
elif [ $1 = "to_nifti" ]; then
    CMD="hb_lawson_to_nifti.py"
elif [ $1 = "bash" ]; then
    CMD="bash"
else
    echo $USAGE
    exit 2
fi

shift

if [ N`uname` = "NDarwin" ]; then
    xhost + $IP
fi
docker run -ti --rm \
    -v `pwd`:/data \
    --env="DISPLAY=$DISP" \
    --volume="/tmp/.X11-unix:/tmp/.X11-unix" \
    --user $(id -u):$(id -g) \
    hb_lawson $CMD $@

if [ N`uname` = "NDarwin" ]; then
    xhost - $IP
fi
