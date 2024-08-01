#!/bin/bash

TRAJ=$1
TPR=$2

if [ -f "OUT_restr.dat" ] ; then
    echo "moved OUT_restr.dat to OUT_restr.dat.bkp"
    mv OUT_restr.dat OUT_restr.dat.bkp
    touch OUT_restr.dat
fi
    touch OUT_restr.dat

echo "protein--ligand" >> OUT_restr.dat
echo " C B A -- a b c" >> OUT_restr.dat

#distance
echo "" >> OUT_restr.dat
echo "" >> OUT_restr.dat
echo "distance" >> OUT_restr.dat
echo 14 | gmx_mpi distance -f $TRAJ -s $TPR -n index_TI.ndx -oall  >> OUT_restr.dat
gmx_mpi analyze -f dist.xvg -dist

#angles
echo "" >> OUT_restr.dat
echo "" >> OUT_restr.dat
echo "angles" >> OUT_restr.dat
echo 15 | gmx_mpi angle -f $TRAJ -n index_TI.ndx -od  >> OUT_restr.dat
mv angdist.xvg angdist_15-BAa.xvg
echo "" >> OUT_restr.dat
echo 16 | gmx_mpi angle -f $TRAJ -n index_TI.ndx -od  >> OUT_restr.dat
mv angdist.xvg angdist_16-Aab.xvg

#dihedrals
echo "" >> OUT_restr.dat
echo "" >> OUT_restr.dat
echo "dihedrals" >> OUT_restr.dat
echo 17 | gmx_mpi angle -f $TRAJ -n index_TI.ndx -od -type dihedral  >> OUT_restr.dat
mv angdist.xvg angdist_17-CBAa.xvg
echo "" >> OUT_restr.dat
echo 18 | gmx_mpi angle -f $TRAJ -n index_TI.ndx -od -type dihedral  >> OUT_restr.dat
mv angdist.xvg angdist_18-Aabc.xvg
echo "" >> OUT_restr.dat
echo 19 | gmx_mpi angle -f $TRAJ -n index_TI.ndx -od -type dihedral  >> OUT_restr.dat
mv angdist.xvg angdist_19-BAab.xvg
