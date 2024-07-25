import pandas as pd 
import numpy as np 
import parmed 
import subprocess
import shlex
import sys
import argparse
import os

from_ff_to_multiego = {
    "OC1": "O1",
    "OC2": "O2",
    "OT1": "O1",
    "OT2": "O2",
}

def parse_arguments():
    parser = argparse.ArgumentParser(description="Process a network topology.")

    # Add the --top flag for the topology path argument
    parser.add_argument(
        "--top",
        metavar="TOPOLOGY_PATH",
        type=str,
        required=True,
        help="Path to the network topology file",
    )

    # Add the --conf flag for the GRO file
    parser.add_argument(
        "--gro",
        metavar="CONF_FILE",
        type=str,
        help="Path to the configuration (GRO) file",
    )

    parser.add_argument(
        "--mego",
        required=False,
        action="store_true",
        help="Path to the configuration (GRO) file",
    )


    # You can add more arguments as needed

    return parser.parse_args()


def get_dihedrals(topology):
    """
    Extracts dihedral angles information from a molecular topology.

    Args:
    - topology (list): List of dihedral atoms information.

    Returns:
    - dihedrals_dataframe (pandas.DataFrame): DataFrame containing dihedral angles data, including atom indices,
      function type, phase, phi_k, and periodicity.
    """
    dihedrals_dataframe = pd.DataFrame(
        {
            "ai": [dihedral.atom1.idx + 1 for dihedral in topology],
            "aj": [dihedral.atom2.idx + 1 for dihedral in topology],
            "ak": [dihedral.atom3.idx + 1 for dihedral in topology],
            "al": [dihedral.atom4.idx + 1 for dihedral in topology],
            #"funct": [dihedral.funct for dihedral in topology],
            #"phase": [dihedral.type.phase for dihedral in topology],
            #"phi_k": [dihedral.type.phi_k for dihedral in topology],
            #"per": [dihedral.type.per for dihedral in topology],
        }
    )
    #dihedrals_dataframe["phi_k"] = dihedrals_dataframe["phi_k"] * 4.184
    return dihedrals_dataframe


def get_impropers(topology):
    """
    Extracts improper torsions information from a molecular topology.

    Args:
    - topology (list): List of improper torsion atoms information.

    Returns:
    - impropers_dataframe (pandas.DataFrame): DataFrame containing improper torsion data, including atom indices,
      function type, psi_eq, and psi_k.
    """
    impropers_dataframe = pd.DataFrame(
        {
            "ai": [improper.atom1.idx + 1 for improper in topology],
            "aj": [improper.atom2.idx + 1 for improper in topology],
            "ak": [improper.atom3.idx + 1 for improper in topology],
            "al": [improper.atom4.idx + 1 for improper in topology],
            "funct": [improper.funct for improper in topology],
            "psi_eq": [improper.type.psi_eq for improper in topology],
            "psi_k": [improper.type.psi_k for improper in topology],
        }
    )
    impropers_dataframe["psi_k"] = impropers_dataframe["psi_k"] * 4.184 * 2
    return impropers_dataframe

def get_pairs(topology):
    """
    Extracts pair information from a molecular topology.

    Args:
    - topology (list): List of pair atoms information.

    Returns:
    - pairs_dataframe (pandas.DataFrame): DataFrame containing pair data, including atom indices, function type, and pair type.
    """
    pairs_dataframe = pd.DataFrame(
        {
            "ai": [pair.atom1.idx + 1 for pair in topology],
            "aj": [pair.atom2.idx + 1 for pair in topology],
            "funct": [pair.funct for pair in topology],
            "type": [pair.type for pair in topology],
        }
    )
    return pairs_dataframe


def get_angles(topology):
    angles_dataframe = pd.DataFrame(
        {
            "ai": [angle.atom1.idx + 1 for angle in topology],
            "aj": [angle.atom2.idx + 1 for angle in topology],
            "ak": [angle.atom3.idx + 1 for angle in topology],
            "funct": [angle.funct for angle in topology],
            "theteq": [angle.type.theteq for angle in topology],
            "k": [angle.type.k for angle in topology],
        }
    )
    angles_dataframe["k"] = angles_dataframe["k"] * 4.184 * 2
    angles_dataframe["k"] = angles_dataframe["k"].map(lambda x: "{:.6e}".format(x))
    return angles_dataframe

###generation of indeces###

def generate_pairs_idx_files(pairs_data_frame, reference_topology,file_content, gro, args):
    dir = "pairs"
    if not os.path.exists(dir):
        try:
            # Create the directory
            os.makedirs(dir)
        except OSError as e:
            print(f"Error creating directory '{dir}': {e}")
    if not os.path.exists(f"{dir}/indeces"):
        try:
            # Create the directory
            os.makedirs(f"{dir}/indeces")
        except OSError as e:
            print(f"Error creating directory : {e}")

    for i in range(len(pairs_data_frame)):
        pairs = [pairs_data_frame["ai"].iloc[i], pairs_data_frame["aj"].iloc[i]]
        atom_names = np.array([reference_topology.atoms[pairs_data_frame["ai"].iloc[i]-1].name, reference_topology.atoms[pairs_data_frame["aj"].iloc[i]-1].name], dtype = str)
        atom_names_mapped = []
        if not args.mego:
            atom_names_mapped = []
            n_H = 0
            for name in atom_names:
                if name[0]=="H":
                    n_H += 1
                else: 
                    if name in list(from_ff_to_multiego):
                        print(f"found {name} to {from_ff_to_multiego[name]}")
                        atom_names_mapped.append(from_ff_to_multiego[name])
                    else: 
                        print("not mapped")
                        atom_names_mapped.append(name)
            if n_H > 0: continue
        entry=f' echo "a {pairs[0]} | a {pairs[1]} \n q " ' 
        p1 = subprocess.Popen(entry,stdout=subprocess.PIPE, shell=True, text=True)
        if args.mego:
            str_out = "_".join(np.sort(np.array(atom_names, dtype = str)))
        else: 
            str_out = "_".join(np.sort(np.array(atom_names_mapped)))
        subprocess.run(shlex.split(f"gmx_mpi make_ndx -f {gro} -o pairs/indeces/index_pair_{str_out}"), stdin=p1.stdout, check=True)



def generate_dihedrals_idx_files(dihedrals_data_frame, reference_topology,file_content, gro, args):
    dir = "dihedrals"
    if not os.path.exists(dir):
        try:
            # Create the directory
            os.makedirs(dir)
        except OSError as e:
            print(f"Error creating directory '{dir}': {e}")
    if not os.path.exists(f"{dir}/indeces"):
        try:
            # Create the directory
            os.makedirs(f"{dir}/indeces")
        except OSError as e:
            print(f"Error creating directory : {e}")

    for i in range(len(dihedrals_data_frame)):
        dihedrals = [dihedrals_data_frame["ai"].iloc[i], dihedrals_data_frame["aj"].iloc[i], dihedrals_data_frame["ak"].iloc[i], dihedrals_data_frame["al"].iloc[i]]
        atom_names = np.array([reference_topology.atoms[dihedrals_data_frame["ai"].iloc[i]-1].name, reference_topology.atoms[dihedrals_data_frame["aj"].iloc[i]-1].name, reference_topology.atoms[dihedrals_data_frame["ak"].iloc[i]-1].name, reference_topology.atoms[dihedrals_data_frame["al"].iloc[i]-1].name], dtype = str)
        atom_names_mapped = []
        if not args.mego:
            atom_names_mapped = []
            n_H = 0
            for name in atom_names:
                if name[0]=="H":
                    n_H += 1
                else: 
                    if name in list(from_ff_to_multiego):
                        print(f"found {name} to {from_ff_to_multiego[name]}")
                        atom_names_mapped.append(from_ff_to_multiego[name])
                    else: 
                        print("not mapped")
                        atom_names_mapped.append(name)
            if n_H > 0: continue
        entry=f' echo "a {dihedrals[0]} | a {dihedrals[1]} | a {dihedrals[2]} | a {dihedrals[3]} \n q " ' 
        p1 = subprocess.Popen(entry,stdout=subprocess.PIPE, shell=True, text=True)
        if args.mego:
            str_out = "_".join(np.sort(np.array(atom_names, dtype = str)))
        else: 
            str_out = "_".join(np.sort(np.array(atom_names_mapped)))
        #subprocess.run(shlex.split(f"gmx_mpi make_ndx -f {gro} -o dihedrals/indeces/index_dihedrals_{str_out}"), stdin=p1.stdout, check=True)
        str_atm_idx_name = "_a_".join(np.array(dihedrals, dtype = str))
        str_atm_idx = "    ".join(np.array(dihedrals, dtype = str))
        modified_content = f"{file_content} \n[a_{str_atm_idx_name}] \n{str_atm_idx}"
        with open(f"dihedrals/indeces/index_dihedrals_{str_out}.ndx", 'w') as output_file:
            output_file.write(modified_content)
    
def generate_impropers_idx_files(impropers_data_frame, reference_topology,file_content ,gro, args):
    dir = "impropers"
    if not os.path.exists(dir):
        try:
            # Create the directory
            os.makedirs(dir)
        except OSError as e:
            print(f"Error creating directory '{dir}': {e}")
    if not os.path.exists(f"{dir}/indeces"):
        try:
            # Create the directory
            os.makedirs(f"{dir}/indeces")
        except OSError as e:
            print(f"Error creating directory : {e}")

    for i in range(len(impropers_data_frame)):
        impropers = [impropers_data_frame["ai"].iloc[i], impropers_data_frame["aj"].iloc[i], impropers_data_frame["ak"].iloc[i], impropers_data_frame["al"].iloc[i]]
        atom_names = np.array([reference_topology.atoms[impropers_data_frame["ai"].iloc[i]-1].name, reference_topology.atoms[impropers_data_frame["aj"].iloc[i]-1].name, reference_topology.atoms[impropers_data_frame["ak"].iloc[i]-1].name, reference_topology.atoms[impropers_data_frame["al"].iloc[i]-1].name], dtype = str)
        if not args.mego:
            atom_names_mapped = []
            n_H = 0
            for name in atom_names:
                if name[0]=="H":
                    n_H += 1
                else: 
                    if name in list(from_ff_to_multiego):
                        print(f"found {name} to {from_ff_to_multiego[name]}")
                        atom_names_mapped.append(from_ff_to_multiego[name])
                    else: 
                        print("not mapped")
                        atom_names_mapped.append(name)
            if n_H > 0: continue
        entry=f' echo "a {impropers[0]} | a {impropers[1]} | a {impropers[2]} | a {impropers[3]} \n q " ' 
        p1 = subprocess.Popen(entry,stdout=subprocess.PIPE, shell=True, text=True)
        if args.mego:
            str_out = "_".join(np.sort(np.array(atom_names, dtype = str)))
        else: 
            str_out = "_".join(np.sort(np.array(atom_names_mapped)))
        #subprocess.run(shlex.split(f"gmx_mpi make_ndx -f {gro} -o impropers/indeces/index_impropers_{str_out}"), stdin=p1.stdout, check=True)
        str_atm_idx_name = "_a_".join(np.array(impropers, dtype = str))
        str_atm_idx = "    ".join(np.array(impropers, dtype = str))
        modified_content = f"{file_content} \n[a_{str_atm_idx_name}] \n{str_atm_idx}"
        with open(f"impropers/indeces/index_impropers_{str_out}.ndx", 'w') as output_file:
            output_file.write(modified_content)

def generate_angle_idx_files(angles_data_frame, reference_topology,file_content,  gro, args):

    dir = "angles"
    if not os.path.exists(dir):
        try:
            # Create the directory
            os.makedirs(dir)
        except OSError as e:
            print(f"Error creating directory '{dir}': {e}")
    if not os.path.exists(f"{dir}/indeces"):
        try:
            # Create the directory
            os.makedirs(f"{dir}/indeces")
        except OSError as e:
            print(f"Error creating directory : {e}")


    for i in range(len(angles_data_frame)):
        angles = [angles_data_frame["ai"].iloc[i], angles_data_frame["aj"].iloc[i], angles_data_frame["ak"].iloc[i]]
        atom_names = np.array([reference_topology.atoms[angles_data_frame["ai"].iloc[i]-1].name, reference_topology.atoms[angles_data_frame["aj"].iloc[i]-1].name, reference_topology.atoms[angles_data_frame["ak"].iloc[i]-1].name], dtype = str)
        
        if not args.mego:
            atom_names_mapped = []
            n_H = 0
            for name in atom_names:
                if name[0]=="H":
                    n_H += 1
                else: 
                    if name in list(from_ff_to_multiego):
                        print(f"found {name} to {from_ff_to_multiego[name]}")
                        atom_names_mapped.append(from_ff_to_multiego[name])
                    else: 
                        print("not mapped")
                        atom_names_mapped.append(name)
            if n_H > 0: continue

        entry=f' echo "a {angles[0]} | a {angles[1]} | a {angles[2]} \n q " ' 

        p1 = subprocess.Popen(entry,stdout=subprocess.PIPE, shell=True, text=True)
        if args.mego:
            str_out = "_".join(np.sort(np.array(atom_names, dtype = str)))
        else: 
            str_out = "_".join(np.sort(np.array(atom_names_mapped)))
        str_atm_idx_name = "_a_".join(np.array(angles, dtype = str))
        str_atm_idx = "    ".join(np.array(angles, dtype = str))
        #subprocess.run(shlex.split(f"gmx_mpi make_ndx -f {gro} -o angles/indeces/index_angle_{str_out}"), stdin=p1.stdout, check=True)
        modified_content = f"{file_content} \n[a_{str_atm_idx_name}] \n{str_atm_idx}"
        with open(f"angles/indeces/index_angle_{str_out}.ndx", 'w') as output_file:
            output_file.write(modified_content)


def main():
    # Parse command-line arguments
    args = parse_arguments()
    
    # Access the topology path
    topology_path = args.top
    reference_topology = parmed.load_file(topology_path)

    #print(reference_topology.atoms[0].name)

    pairs_data_frame = get_pairs(reference_topology.adjusts)
    angles_data_frame = get_angles(reference_topology.angles)
    dihedrals_data_frame = get_dihedrals(reference_topology.dihedrals)
    impropers_data_frame = get_impropers(reference_topology.impropers)
    with open("index.ndx", 'r') as file:
        file_content = file.read()
    generate_pairs_idx_files(pairs_data_frame, reference_topology,file_content,  args.gro, args)
    generate_angle_idx_files(angles_data_frame, reference_topology,file_content,  args.gro, args)
    generate_dihedrals_idx_files(dihedrals_data_frame, reference_topology,file_content, args.gro, args)
    generate_impropers_idx_files(impropers_data_frame, reference_topology,file_content, args.gro, args)

if __name__ == "__main__":
    main()

