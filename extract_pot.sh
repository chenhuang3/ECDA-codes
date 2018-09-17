#!/bin/bash

module purge 
module load intel-openmpi

###cut3d=/home/chuang3/work_fsu/codes/abinit-7.10.5/src/98_main/cut3d
cut3d=/gpfs/research/huang/chenh/abinit-7.10.5/src/98_main/cut3d

echo "go_POT
   1
   1
   5
   up.dat
0  " | $cut3d


echo "go_POT
   1
   2
   5
   down.dat
0  " | $cut3d


ls -lh  up.dat
ls -lh  down.dat

cat up.dat > ../init_guess_total_vks.dat
cat down.dat >> ../init_guess_total_vks.dat

ls -lh ../init_guess_total_vks.dat

echo ""
echo ""
echo "done!"
echo ""
echo ""
echo "DONOT FORGET TO COPY VPSP.Dat for xcpatch as well"
echo ""
echo ""
echo ""
