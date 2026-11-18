#!/bin/bash

start=10
end=1530
step=20
wd=`pwd`
coords=(98 136 102 118 154 168 172 174 166 170 178 220 224 232 234 242 244 354 356 368 374 386 376 370 388 462 496 476 478 500 504 506 508 510 514 502 528)
# Loop through the values and create directories

for crdnt in $(seq $start $step $end); do
    echo "$crdnt"
    rm -r "coord_$crdnt"
    mkdir -p "coord_$crdnt"
    cd "coord_$crdnt" || exit
    cp ../run.sh .
    sbatch --export=crdnt=$crdnt,wd=$wd run.sh
    cd "$wd" || exit
done
