import pandas as pd 
import numpy as np 
import subprocess
import shlex
import sys
import argparse
import os
import argparse
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser(description="Process a network topology.")

parser.add_argument(
    "--mego",
    metavar="mego path",
    type=str,
    required=True,
    help="Path to the network topology file",
)
parser.add_argument(
    "--aa",
    metavar="all_atom path",
    type=str,
    required=True,
    help="Path to the network topology file",
)

parser.add_argument(
    "--comparison",
    metavar="file_for_comparison path",
    type=str,
    required=True,
    help="Path to the network topology file",
)
args = parser.parse_args()

dir = f"{args.comparison}"
if not os.path.exists(f"{dir}"):
    try:
        # Create the directory
        os.makedirs(dir)
    except OSError as e:
        print(f"Error creating directory '{dir}': {e}")

#Pairs
def plot_pairs(args):
    print("starting Pairs")
    dir = f"{args.comparison}/pairs"
    if not os.path.exists(f"{dir}"):
        try:
            # Create the directory
            os.makedirs(dir)
        except OSError as e:
            print(f"Error creating directory '{dir}': {e}")

    distr_mego = [x for x in os.listdir(f"{args.mego}/pairs/pairs_distr") if "#" not in x]
    distr_aa   = [x for x in os.listdir(f"{args.aa}/pairs/pairs_distr") if "#" not in x]
    common_elements = list(set(distr_mego) & set(distr_aa))

    for i in range(len(common_elements)):
        name = common_elements[i].split("pairs_d_")[1].split(".xvg")[0]
        angles_mego, dis_mego = np.loadtxt(f"{args.mego}/pairs/pairs_distr/{common_elements[i]}", unpack=True, comments=["#", "@"])
        angles_aa, dis_aa = np.loadtxt(f"{args.aa}/pairs/pairs_distr/{common_elements[i]}", unpack=True, comments=["#", "@"])
        fig, ax = plt.subplots(1,1)
        bins=30
        ax.hist(dis_mego, label="mego", density=True, alpha = 0.6, bins=bins)
        ax.hist(dis_aa, label="train", density=True, alpha = 0.6, bins=bins)
        ax.legend()
        fig.savefig(f"{dir}/{name}.png")

#Angles
def plot_angles(args):
    print("starting Angles")

    dir = f"{args.comparison}/angles"
    if not os.path.exists(f"{dir}"):
        try:
            # Create the directory
            os.makedirs(dir)
        except OSError as e:
            print(f"Error creating directory '{dir}': {e}")

    distr_mego = [x for x in os.listdir(f"{args.mego}/angles/angles_distr") if "#" not in x]
    distr_aa   = [x for x in os.listdir(f"{args.aa}/angles/angles_distr") if "#" not in x]
    common_elements = list(set(distr_mego) & set(distr_aa))

    for i in range(len(common_elements)):
        name = common_elements[i].split("ang_dis_")[1].split(".xvg")[0]
        angles_mego, dis_mego = np.loadtxt(f"{args.mego}/angles/angles_distr/{common_elements[i]}", unpack=True, comments=["#", "@"])
        angles_aa, dis_aa = np.loadtxt(f"{args.aa}/angles/angles_distr/{common_elements[i]}", unpack=True, comments=["#", "@"])

        fig, ax = plt.subplots(1,1)
        ax.plot(angles_mego, dis_mego, label="mego")
        ax.plot(angles_aa, dis_aa, label="train")
        ax.grid()
        ax.legend()
        fig.savefig(f"{dir}/{name}.png")

#Dihedrals
def plot_dihedrals(args):
    print("starting Dihedrals")

    dir = f"{args.comparison}/dihedrals"
    if not os.path.exists(f"{dir}"):
        try:
            # Create the directory
            os.makedirs(dir)
        except OSError as e:
            print(f"Error creating directory '{dir}': {e}")

    distr_mego = [x for x in os.listdir(f"{args.mego}/dihedrals/dihedrals_distr/") if "#" not in x]
    distr_aa   = [x for x in os.listdir(f"{args.aa}/dihedrals/dihedrals_distr/") if "#" not in x]

    common_elements = list(set(distr_mego) & set(distr_aa))
    for i in range(len(common_elements)):
        name = common_elements[i].split("dih_dis")[1].split(".xvg")[0]
        angles_mego, dis_mego = np.loadtxt(f"{args.mego}/dihedrals/dihedrals_distr/{common_elements[i]}", unpack=True, comments=["#", "@"])
        angles_aa, dis_aa = np.loadtxt(f"{args.aa}/dihedrals/dihedrals_distr/{common_elements[i]}", unpack=True, comments=["#", "@"])

        fig, ax = plt.subplots(1,1)
        ax.plot(angles_mego, dis_mego, label="mego")
        ax.plot(angles_aa, dis_aa, label="train")
        ax.legend()
        ax.grid()
        fig.savefig(f"{dir}/{name}.png")

#Impropers
def plot_impropers(args):
    print("starting Impropers")

    dir = f"{args.comparison}/impropers"
    if not os.path.exists(f"{dir}"):
        try:
            # Create the directory
            os.makedirs(dir)
        except OSError as e:
            print(f"Error creating directory '{dir}': {e}")

    distr_mego = [x for x in os.listdir(f"{args.mego}/impropers/impropers_distr/") if "#" not in x]
    distr_aa   = [x for x in os.listdir(f"{args.aa}/impropers/impropers_distr/") if "#" not in x]

    common_elements = list(set(distr_mego) & set(distr_aa))

    for i in range(len(common_elements)):
        name = common_elements[i].split("impropers_dis")[1].split(".xvg")[0]
        angles_mego, dis_mego = np.loadtxt(f"{args.mego}/impropers/impropers_distr/{common_elements[i]}", unpack=True, comments=["#", "@"])
        angles_aa, dis_aa = np.loadtxt(f"{args.aa}/impropers/impropers_distr/{common_elements[i]}", unpack=True, comments=["#", "@"])
        fig, ax = plt.subplots(1,1)
        ax.plot(angles_mego, dis_mego, label="mego")
        ax.plot(angles_aa, dis_aa, label="train")
        ax.legend()
        ax.grid()
        fig.savefig(f"{dir}/{name}.png")


plot_angles(args)
plot_dihedrals(args)
plot_impropers(args)
plot_pairs(args)
