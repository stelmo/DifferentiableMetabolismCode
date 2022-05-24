#!/bin/bash

# master_ids=("pgi8#B2" "pgi1#B2" "pts3#B1" "pgi1#B1" "pgi5#B1" "tpi1#B1" "sdh1#B1" "tpi2#B1" "pgi4#B2" "tpi2#B2" "tpi4#B2" "pgi2#B2" "sdh1#B2" "sdh3#B2" "tpi3#B2" "WT2#B1" "pgi7#B2" "WT2#B2" "pts2#B2" "pgi4#B1" "tpi3#B1" "pts3#B2" "WT1#B2" "pgi7#B1" "pts4#B1" "pts4#B2" "tpi1#B2" "pgi5#B2" "pgi6#B1" "pts1#B2" "pts2#B1" "pgi3#B1" "tpi4#B1" "pts1#B1" "pgi2#B1" "sdh2#B1" "WT1#B1"  "pgi6#B2" "pgi3#B2" "pgi8#B1" "sdh3#B1" "sdh2#B2")

for master_id in "${master_ids[@]}"
do
    echo launch $master_id
    `nohup julia ./analyses/gd_gecko/cluster.jl "$master_id" >  $master_id.out &`
done
