import pandas as pd 
import numpy as np 
import parmed 
import subprocess
import shlex
import sys
import argparse
import os
import argparse


parser = argparse.ArgumentParser(description="Process a network topology.")

parser.add_argument(
    "--xtc",
    metavar="trajectory path",
    type=str,
    required=True,
    help="Path to the network topology file",
)
parser.add_argument("-d", "--dihedrals", action="store_true", help="Include dihedral entries")
parser.add_argument("-p", "--pairs", action="store_true", help="Include pairs entries")
parser.add_argument("-a", "--angles", action="store_true", help="Include angles entries")
parser.add_argument("-i", "--impropers", action="store_true", help="Include impropers entries")

args = parser.parse_args()
#if none of them is specified calculate all of them
if args.dihedrals == False and args.pairs == False and args.angles == False and args.impropers == False:
    args.dihedrals = True
    args.pairs = True
    args.angles = True
    args.impropers = True

#DISTANCES
def calc_pairs():
    idxs = os.listdir("pairs/indeces")
    dir = "pairs/pairs_distr"
    if not os.path.exists(f"{dir}"):
        try:
            # Create the directory
            os.makedirs(dir)
        except OSError as e:
            print(f"Error creating directory '{dir}': {e}")

    for idx in idxs:
        name = idx.split("index_pair_")[1].split(".ndx")[0]

        entry=f' echo "3" ' 
        print(entry)
        p1 = subprocess.Popen(entry,stdout=subprocess.PIPE, shell=True, text=True)
        subprocess.run(shlex.split(f"gmx_mpi distance -f {args.xtc} -n pairs/indeces/{idx} -oall {dir}/pairs_d_{name}"), stdin=p1.stdout, check=True)


#ANGLES
def calc_angles():

    idxs = os.listdir("angles/indeces")
    dir = "angles/angles_distr"
    if not os.path.exists(f"{dir}"):
        try:
            # Create the directory
            os.makedirs(dir)
        except OSError as e:
            print(f"Error creating directory '{dir}': {e}")

    for idx in idxs:
        name = idx.split("index_angle_")[1].split(".ndx")[0]

        entry=f' echo "3" ' 
        print(entry)
        p1 = subprocess.Popen(entry,stdout=subprocess.PIPE, shell=True, text=True)
        subprocess.run(shlex.split(f"gmx_mpi angle -f {args.xtc} -n angles/indeces/{idx} -binwidth 1 -od {dir}/ang_dis_{name}"), stdin=p1.stdout, check=True)

#Dihedrals
def calc_dihedrals():

    idxs = os.listdir("dihedrals/indeces")
    dir = "dihedrals/dihedrals_distr"
    if not os.path.exists(f"{dir}"):
        try:
            # Create the directory
            os.makedirs(dir)
        except OSError as e:
            print(f"Error creating directory '{dir}': {e}")
    for idx in idxs:
        name = idx.split("index_dihedrals_")[1].split(".ndx")[0]

        entry=f' echo "3" ' 
        print(entry)
        p1 = subprocess.Popen(entry,stdout=subprocess.PIPE, shell=True, text=True)
        subprocess.run(shlex.split(f"gmx_mpi angle -f {args.xtc} -n dihedrals/indeces/{idx} -binwidth 3 -od {dir}/dih_dis_{name} --type dihedral"), stdin=p1.stdout, check=True)

#Impropers
def calc_impropers():
    idxs = os.listdir("impropers/indeces")
    dir = "impropers/impropers_distr"
    if not os.path.exists(f"{dir}"):
        try:
            # Create the directory
            os.makedirs(dir)
        except OSError as e:
            print(f"Error creating directory '{dir}': {e}")
    for idx in idxs:
        name = idx.split("index_impropers_")[1].split(".ndx")[0]

        entry=f' echo "3" ' 
        print(entry)
        p1 = subprocess.Popen(entry,stdout=subprocess.PIPE, shell=True, text=True)
        subprocess.run(shlex.split(f"gmx_mpi angle -f {args.xtc} -n impropers/indeces/{idx} -binwidth 5 -od {dir}/impropers_dis_{name} --type dihedral"), stdin=p1.stdout, check=True)

#calculate the specified fields
if args.dihedrals: calc_dihedrals()
if args.angles:    calc_angles()
if args.impropers: calc_impropers()
if args.pairs:     calc_pairs()
