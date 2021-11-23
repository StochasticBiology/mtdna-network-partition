# compile simulation code
gcc -o3 network-sim.c -lm -o network-sim.ce

# run different parameterisations in parallel
./network-sim.ce 0.1 4 > tmp1 &
./network-sim.ce 0.1 16 > tmp2 &
./network-sim.ce 0.1 64 > tmp3 &
./network-sim.ce 0.5 4 > tmp4 &
./network-sim.ce 0.5 16 > tmp5 &
./network-sim.ce 0.5 64 > tmp6 &
./network-sim.ce 0.9 4 > tmp7 &
./network-sim.ce 0.9 16 > tmp8 &
./network-sim.ce 0.9 64 > tmp9 &

# after these jobs have terminated, run analyse.sh
# ./analyse.sh
