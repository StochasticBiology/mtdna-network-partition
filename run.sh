# compile simulation code
gcc -o3 network-sim.c -lm -o network-sim.ce

# run different parameterisations in parallel
./network-sim.ce 0.25 0.5 4 > tmp1 &
./network-sim.ce 0.25 0.5 16 > tmp2 &
./network-sim.ce 0.25 0.5 64 > tmp3 &
./network-sim.ce 0.1 0.5 4 > tmp4 &
./network-sim.ce 0.1 0.5 16 > tmp5 &
./network-sim.ce 0.1 0.5 64 > tmp6 &
./network-sim.ce 0. 0.5 4 > tmp7 &
./network-sim.ce 0. 0.5 16 > tmp8 &
./network-sim.ce 0. 0.5 64 > tmp9 &
./network-sim.ce 0.8 0.5 4 > tmp10 &
./network-sim.ce 0.8 0.5 16 > tmp11 &
./network-sim.ce 0.8 0.5 64 > tmp12 &

# after these jobs have terminated, run analyse.sh
# ./analyse.sh

./network-sim.ce 0.25 0.1 4 > tmp1 &
./network-sim.ce 0.25 0.1 16 > tmp2 &
./network-sim.ce 0.25 0.1 64 > tmp3 &
./network-sim.ce 0.1 0.1 4 > tmp4 &
./network-sim.ce 0.1 0.1 16 > tmp5 &
./network-sim.ce 0.1 0.1 64 > tmp6 &
./network-sim.ce 0. 0.1 4 > tmp7 &
./network-sim.ce 0. 0.1 16 > tmp8 &
./network-sim.ce 0. 0.1 64 > tmp9 &
./network-sim.ce 0.8 0.1 4 > tmp10 &
./network-sim.ce 0.8 0.1 16 > tmp11 &
./network-sim.ce 0.8 0.1 64 > tmp12 &
