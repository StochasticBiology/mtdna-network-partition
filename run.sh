#!/bin/bash
# compile simulation code
gcc -o3 network-sim.c -lm -o network-sim.ce

./network-sim.ce --snapshots 0.5 0.5 4 0 0 0 0 > tmps1 &
./network-sim.ce --snapshots 0.5 0.5 4 0.5 0.5 0 0 > tmps2 &
./network-sim.ce --snapshots 0.5 0.5 4 1 0 0 0 > tmps3 &
./network-sim.ce --snapshots 0.5 0.5 16 0.5 0.5 0 0 > tmps4 &
./network-sim.ce --snapshots 0.5 0.5 64 0.5 0.5 0 0 > tmps5 &
./network-sim.ce --snapshots 0.5 0.5 16 1 1 0 0 > tmps6 &
./network-sim.ce --snapshots 0.5 0.5 16 1 1 0.04 0 > tmps7 &
./network-sim.ce --snapshots 0.5 0.5 16 1 1 0.1 0 > tmps8 &
./network-sim.ce --snapshots 0.5 0.5 16 1 1 0. 0.1 > tmps9 &

# prop = {0.10,0.25,0.50} for 10%, 25%, 50%, respectively

# run different parameterisations in parallel
# h = 0.5 for different prop, seeds
./network-sim.ce --simulate 0.1 0.5 4 > tmp1 &
./network-sim.ce --simulate 0.1 0.5 16 > tmp2 &
./network-sim.ce --simulate 0.1 0.5 64 > tmp3 &
./network-sim.ce --simulate 0.25 0.5 4 > tmp4 &
./network-sim.ce --simulate 0.25 0.5 16 > tmp5 &
./network-sim.ce --simulate 0.25 0.5 64 > tmp6 &
./network-sim.ce --simulate 0.5 0.5 4 > tmp7 &
./network-sim.ce --simulate 0.5 0.5 16 > tmp8 &
./network-sim.ce --simulate 0.5 0.5 64 > tmp9 &

# h = 0.1 for different prop, seeds
./network-sim.ce --simulate 0.1 0.1 4 > tmp1 &
./network-sim.ce --simulate 0.1 0.1 16 > tmp2 &
./network-sim.ce --simulate 0.1 0.1 64 > tmp3 &
./network-sim.ce --simulate 0.25 0.1 4 > tmp4 &
./network-sim.ce --simulate 0.25 0.1 16 > tmp5 &
./network-sim.ce --simulate 0.25 0.1 64 > tmp6 &
./network-sim.ce --simulate 0.5 0.1 4 > tmp7 &
./network-sim.ce --simulate 0.5 0.1 16 > tmp8 &
./network-sim.ce --simulate 0.5 0.1 64 > tmp9 &

# after these jobs have terminated, run analyse.sh
#./analyse-all.sh
