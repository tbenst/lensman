#!/usr/bin/env bash
if [[ $# -eq 0 ]] ; then
    echo 'Usage: transfer.sh REMOTE_FISHDIR REMOTE_TIFDIR'
    echo 'e.g. scripts/transfer.sh 2021-01-19_chrmine_kv2.1_h2b6s_7dpf/fish1_chrmine TSeries-1024cell-32concurrent-4power-046'
    exit 0
fi

fishDir=$1
tifDir=$2
mkdir -p /data/dlab/b115/$fishDir/$tifDir

rsync -avP \
    tbenst@sher:"/oak/stanford/groups/deissero/users/tyler/b115/$fishDir/$tifDir.ty.h5" \
    /data/dlab/b115/$fishDir/

# copy voltage recording
scp \
    tbenst@sher:"/oak/stanford/groups/deissero/users/tyler/b115/$fishDir/$tifDir/${tifDir}_Cycle00001_VoltageRecording*.csv" \
    /data/dlab/b115/$fishDir/$tifDir/

# copy prairie view settings
scp \
    tbenst@sher:"/oak/stanford/groups/deissero/users/tyler/b115/$fishDir/$tifDir/${tifDir}.xml" \
    /data/dlab/b115/$fishDir/$tifDir/

# copy SLM stim files
scp \
    tbenst@sher:"/oak/stanford/groups/deissero/users/tyler/b115/$fishDir/*.mat" \
    /data/dlab/b115/$fishDir/
scp \
    tbenst@sher:"/oak/stanford/groups/deissero/users/tyler/b115/$fishDir/*.txt" \
    /data/dlab/b115/$fishDir/

echo "Done! see /data/dlab/b115/$fishDir/$tifDir/"