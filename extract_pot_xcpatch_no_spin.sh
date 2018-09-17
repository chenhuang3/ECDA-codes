#!/usr/bin/bash

module purge 
module load intel-openmpi

cut3d="mpirun -np 1   /panfs/storage.local/scs/huang/chenh/abinit-7.10.5/src/98_main/cut3d"

echo "go_POT
   1
   5
   up.dat
0  " | $cut3d

ls -lh  up.dat

cat up.dat    > ../init_guess_total_vks.dat

######awk '{getline t<"hxc_up.dat";  print $0-t}' up.dat > ../global_vpsp.dat

ls -lh ../init_guess_total_vks.dat

echo ""
echo ""
echo ""
echo "done!"
echo ""
echo ""
echo "cp ../init_guess_total_vks.dat  ../init_guess_total_vks.dat to working folders"
echo ""
echo "Get vpsp.dat from the global_system, by running and stopping the xcpatch, " 
echo "once you run global_system, just copy the vpsp.dat to this folder as global_vpsp.dat"
echo ""
echo ""
echo ""

rm up.dat
rm down.dat
