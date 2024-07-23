import numpy as np
import os
import math
import json
import argparse


def read_json_file(filepath):
    with open(filepath, 'r') as file:
        data = json.load(file)
    atom_idx = data['atom_idx']
    x_0 = data['x_0']
    k = data['k']
    return atom_idx, x_0, k

def calc_restr(r0, theta0, theta1, K_r, K_t, K_t2, K_phi, K_phi_2, K_phi_3, T):
    theta0 = math.radians(theta0)  # convert angle from degrees to radians --> math.sin() wants radians
    theta1 = math.radians(theta1)  # convert angle from degrees to radians --> math.sin() wants radians
    arg =(
        (8.0 * math.pi**2.0 * V) / (r0**2.0 * math.sin(theta0) * math.sin(theta1))   *   (( ( K_r*K_t* K_phi*K_t2*K_phi_2*K_phi_3 )**0.5 ) / ( (2.0 * math.pi * K * T)**(3) ) )
    )

    dG = - K * T * math.log(arg)
    dG_sim = 0#K*T*np.log(1) no simmetry in the conformation of the molecule that can be broken by the restraint

    return np.abs(dG)+np.abs(dG_sim)

K = 8.314472*0.001  # Gas constant in kJ/mol/K
V = 1.66            # standard volume in nm^3

def main():
    parser = argparse.ArgumentParser(description='Generate formatted outputs for molecular interaction setups.')
    parser.add_argument('--json', type=str, required=True, help='json file with data')
    parser.add_argument('--T', type=float, required=True, help='temperature of the simulation')
    args = parser.parse_args()

    input_numbers, x0, k = read_json_file(args.json)

    dirs = [i for i in os.listdir() if "run" in i]

    dgs = []
    for d in dirs:
        if os.path.exists(f"{d}/all_md"):
            print(f"all_md in {d}")
            if os.path.exists(f"{d}/all_md/barint.xvg"):
                print("reading barint.xvg")
                l, dg = np.loadtxt(f"{d}/all_md/barint.xvg", unpack=True, comments=["#", "@"])
                dg_r = calc_restr(x0[0], x0[1] , x0[2], k[0], k[1], k[2], k[3], k[4], k[5], args.T)
                dgs.append(dg[-1]*K*args.T - dg_r)
            else: 
                print("""
File barint.xvg not found. 
calculate with 
gmx_mpi bar -f md* -o -oi
                      """)

    dgs = np.array(dgs)
    print(f"DGs: {dgs}")
    print(f"Average Free energy: {np.mean(dgs)}")
    print(f"Free energy std: {np.std(dgs)}")
    print(f"Std of the mean: {np.std(dgs)/np.sqrt(len(dgs)-1)}")
    print()
    print(f"<DG> = {np.mean(dgs)} +- {np.std(dgs)/np.sqrt(len(dgs)-1)}")

    with open("DG.dat", 'w') as file:
        # Write content to the file
        file.write(f"DGs: {dgs}")
        file.write(f"<DG> = {np.mean(dgs)} +- {np.std(dgs)/np.sqrt(len(dgs)-1)}")

if __name__ == "__main__":
    main()