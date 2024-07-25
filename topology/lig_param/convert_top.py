import parmed as pmd 
import datetime
import pandas as pd 
import numpy as np
import argparse
from parmed import unit as u
import os
import warnings

warnings.filterwarnings("ignore")
date = datetime.datetime.now()
def make_header(args):
    header = f"""
;
;  Modified topology: removed hydrogen and reindexed all indeces accordingly
;  Starting topology input system: input/{args.system}
;  Ending topology output: input{args.system} 
;  At date: {date.strftime("%c")}
;  Scale masses: {args.scale_masses}
;  Scale bonds : {args.scale_bonds}

"""
    return header

from_ff_to_multiego = {
    "OC1": "O1",
    "OC2": "O2",
    "OT1": "O1",
    "OT2": "O2",
}

def get_atoms_types(topology):

    atoms_dataframe = pd.DataFrame(
        {
            "name": [atom.name for atom in topology],
            "type": [atom.type for atom in topology],
            "at.num": [atom.atomic_number for atom in topology],
            "mass": [atom.mass for atom in topology],
            "charge": [0.0 for atom in topology],
            "ptype" : ["A" for atom in topology],
            "atom": [atom.idx + 1 for atom in topology],
            "sigma": [atom.sigma*0.1 for atom in topology],
            "epsilon": [atom.epsilon*4.184 for atom in topology],
        }
    )
    return atoms_dataframe

def get_atoms(topology):

    atoms_dataframe = pd.DataFrame(
        {
            "nr": [atom.idx + 1 for atom in topology],
            "type": [atom.type for atom in topology],
            "resnr": [atom.residue.idx + 1 for atom in topology],
            "residue": [atom.residue.name for atom in topology],
            "name": [atom.name for atom in topology],
            "atom": [atom.idx + 1 for atom in topology],
            "charge": [0.0 for atom in topology],
            "mass": [atom.mass for atom in topology],
        }
    )
    return atoms_dataframe

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
            "funct": [dihedral.funct for dihedral in topology],
            "phase": [ [ multi_dih_type.phase for multi_dih_type in dihedral.type ] if topology[0].funct==9 else dihedral.type.phase for dihedral in topology],
            "phi_k": [ [ multi_dih_type.phi_k * 4.184 for multi_dih_type in dihedral.type ] if topology[0].funct==9 else dihedral.type.phi_k * 4.184 for dihedral in topology],
            "per": [ [ multi_dih_type.per for multi_dih_type in dihedral.type ] if topology[0].funct==9 else dihedral.type.per for dihedral in topology],
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

    if topology[0].funct == 1:
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
    
    if topology[0].funct == 5:
        angles_dataframe = pd.DataFrame(
            {
                "ai": [angle.atom1.idx + 1 for angle in topology],
                "aj": [angle.atom2.idx + 1 for angle in topology],
                "ak": [angle.atom3.idx + 1 for angle in topology],
                "funct": [angle.funct for angle in topology],
                "theteq": [angle.type.theteq for angle in topology],
                "k": [angle.type.k for angle in topology],
                #"s0": [angle.type.k for angle in topology],
                #"Kub": [angle.type.k for angle in topology],
            }
        )
        angles_dataframe["k"] = angles_dataframe["k"] * 4.184 * 2 * 10
        angles_dataframe["k"] = angles_dataframe["k"].map(lambda x: "{:.6e}".format(x))
    
    return angles_dataframe

def get_bonds(topology):

    """
    Extracts bonds  information from a molecular topology.

    Args:
    - topology (list): List of dihedral atoms information.

    Returns:
    - bonds_dataframe (pandas.DataFrame): DataFrame containing dihedral angles data, including atom indices,
      function type, phase, phi_k, and periodicity.
    """
    bonds_dataframe = pd.DataFrame(
        {
            "ai": [bond.atom1.idx + 1 for bond in topology],
            "aj": [bond.atom2.idx + 1 for bond in topology],
            "funct": [bond.funct for bond in topology],
            "c0": [bond.type.req for bond in topology],
            "c1": [bond.type.k for bond in topology],
        }
    )
    bonds_dataframe["c1"] = bonds_dataframe["c1"] * 4.184 *100 * 2 # reconvert to unit of measurment of des amber
    #bonds_dataframe["c1"] = bonds_dataframe["c1"] * 4.184/100     # convert to kJ/mol/nm^2
    bonds_dataframe["c0"] = bonds_dataframe["c0"] * 0.1 #reconvert to nm
    return bonds_dataframe

def convert_bonds(bonds, h_atom_num, mapping_dict):

    bonds_noH = bonds.copy()
    mask1 = bonds_noH['ai'].isin(h_atom_num)
    mask2 = bonds_noH['aj'].isin(h_atom_num)
    mask = np.logical_or(mask1, mask2)

    #remove Hydrogens
    bonds_noH = bonds_noH[~mask]

    #map indeces
    mapped_indices1 = [mapping_dict[index] for index in bonds_noH["ai"]]
    mapped_indices2 = [mapping_dict[index] for index in bonds_noH["aj"]]

    #reassign indeces to bonds
    bonds_noH["ai"] = mapped_indices1
    bonds_noH["aj"] = mapped_indices2

    return bonds_noH

def convert_pairs(pairs, h_atom_num, mapping_dict):

    pairs_noH = pairs.copy()
    mask1 = pairs_noH['ai'].isin(h_atom_num)
    mask2 = pairs_noH['aj'].isin(h_atom_num)
    mask = np.logical_or(mask1, mask2)

    #remove Hydrogens
    pairs_noH = pairs_noH[~mask]

    #map indeces
    mapped_indices1 = [mapping_dict[index] for index in pairs_noH["ai"]]
    mapped_indices2 = [mapping_dict[index] for index in pairs_noH["aj"]]

    #reassign indeces to bonds
    pairs_noH["ai"] = mapped_indices1
    pairs_noH["aj"] = mapped_indices2

    return pairs_noH

def convert_angles(angles, h_atom_num, mapping_dict):

    angles_noH = angles.copy()
    mask1 = angles_noH['ai'].isin(h_atom_num)
    mask2 = angles_noH['aj'].isin(h_atom_num)
    mask3 = angles_noH['ak'].isin(h_atom_num)
    mask = np.logical_or(np.logical_or(mask1, mask2), mask3)

    #remove Hydrogens
    angles_noH = angles_noH[~mask]

    #map indeces
    mapped_indices1 = [mapping_dict[index] for index in angles_noH["ai"]]
    mapped_indices2 = [mapping_dict[index] for index in angles_noH["aj"]]
    mapped_indices3 = [mapping_dict[index] for index in angles_noH["ak"]]

    #reassign indeces to bonds
    angles_noH["ai"] = mapped_indices1
    angles_noH["aj"] = mapped_indices2
    angles_noH["ak"] = mapped_indices3

    return angles_noH

def convert_dihedrals(dihedrals, h_atom_num, mapping_dict):

    dihedrals_noH = dihedrals.copy()
    mask1 = dihedrals_noH['ai'].isin(h_atom_num)
    mask2 = dihedrals_noH['aj'].isin(h_atom_num)
    mask3 = dihedrals_noH['ak'].isin(h_atom_num)
    mask4 = dihedrals_noH['al'].isin(h_atom_num)

    mask = np.logical_or(np.logical_or(np.logical_or(mask1, mask2), mask3), mask4)

    #remove Hydrogens
    dihedrals_noH = dihedrals_noH[~mask]

    #map indeces
    mapped_indices1 = [mapping_dict[index] for index in dihedrals_noH["ai"]]
    mapped_indices2 = [mapping_dict[index] for index in dihedrals_noH["aj"]]
    mapped_indices3 = [mapping_dict[index] for index in dihedrals_noH["ak"]]
    mapped_indices4 = [mapping_dict[index] for index in dihedrals_noH["al"]]

    #reassign indeces to bonds
    dihedrals_noH["ai"] = mapped_indices1
    dihedrals_noH["aj"] = mapped_indices2
    dihedrals_noH["ak"] = mapped_indices3
    dihedrals_noH["al"] = mapped_indices4

    if dihedrals_noH["funct"].to_numpy()[0]==9:
        cols = dihedrals_noH.columns
        dih_appo = pd.DataFrame(columns=cols)
        for i in range(len(dihedrals_noH)):
            d = dihedrals_noH.iloc[i].copy()
            for i_appo in range(len(d["phase"])):
                i_dih = pd.DataFrame({
                    "ai"    : [d["ai"]],
                    "aj"    : [d["aj"]],
                    "ak"    : [d["ak"]],
                    "al"    : [d["al"]],
                    "funct" : [1],
                    "phase" : [d["phase"][i_appo]],
                    "phi_k" : [d["phi_k"][i_appo]],
                    "per"   : [d["per"][i_appo]],
                    }
                )

                dih_appo = pd.concat([dih_appo, i_dih])
        dihedrals_noH = dih_appo.copy()

    return dihedrals_noH

def convert_impropers(dihedrals, h_atom_num, mapping_dict):

    dihedrals_noH = dihedrals.copy()
    mask1 = dihedrals_noH['ai'].isin(h_atom_num)
    mask2 = dihedrals_noH['aj'].isin(h_atom_num)
    mask3 = dihedrals_noH['ak'].isin(h_atom_num)
    mask4 = dihedrals_noH['al'].isin(h_atom_num)

    mask = np.logical_or(np.logical_or(np.logical_or(mask1, mask2), mask3), mask4)

    #remove Hydrogens
    dihedrals_noH = dihedrals_noH[~mask]

    #map indeces
    mapped_indices1 = [mapping_dict[index] for index in dihedrals_noH["ai"]]
    mapped_indices2 = [mapping_dict[index] for index in dihedrals_noH["aj"]]
    mapped_indices3 = [mapping_dict[index] for index in dihedrals_noH["ak"]]
    mapped_indices4 = [mapping_dict[index] for index in dihedrals_noH["al"]]

    #reassign indeces to bonds
    dihedrals_noH["ai"] = mapped_indices1
    dihedrals_noH["aj"] = mapped_indices2
    dihedrals_noH["ak"] = mapped_indices3
    dihedrals_noH["al"] = mapped_indices4

    return dihedrals_noH

def get_H_mapping(atoms, args):

    #copy the atom type and remove hydrogens
    atoms_noH = atoms.copy()

    #print(len(atoms_noH["mass"]))
    atoms_noH = atoms_noH[~atoms_noH['name'].str.startswith('H')]
    atoms_H   = atoms[atoms['name'].str.startswith('H')]

    #Calculate new masses: redistribute H mass on heavy atom
    if args.scale_masses: 
        mass_H_removed = np.zeros(len(atoms_noH["mass"]))
        mass_appo,idx_heavy_atom,idx_appo = 0,0,0
        for i in range(len(atoms)):
            if atoms["name"].to_numpy()[i][0]!='H':
                if i > 0 :
                    mass_H_removed[idx_appo] = atoms["mass"].to_numpy()[idx_heavy_atom] + mass_appo
                    idx_appo += 1

                idx_heavy_atom = i
                mass_appo = 0
            else:
                mass_appo += atoms["mass"].to_numpy()[i]
        if idx_appo < len(mass_H_removed): mass_H_removed[idx_appo] = atoms["mass"].to_numpy()[idx_heavy_atom] + mass_appo

        atoms_noH["mass"] = mass_H_removed

    #save the indexes of the hydrogens
    h_atom_num = np.array(atoms_H.index) + 1

    #save the old numbering
    old_enumeration = np.array(atoms_noH.index) + 1
    #reset indexes
    atoms_noH = atoms_noH.reset_index(drop=True)
    #save new numbering (without hydrogens)
    new_enumeration = np.array(atoms_noH.index) + 1

    #create the mapping function to map the "with H enumeration" to the "without H enumeration"
    mapping_dict = dict(zip(old_enumeration, new_enumeration))

    return mapping_dict, h_atom_num, atoms_noH

def dataframe_to_write(df, header=True, col_space = []):
    """
    Returns a stringified and formated dataframe and a message if the dataframe is empty.

    Parameters
    ----------
    df : pd.DataFrame
        The input dataframe

    Returns
    -------
    The stringified dataframe
    """
    if df.empty:
        # TODO insert and improve the following warning
        print("A topology parameter is empty. Check the reference topology.")
        return "; The following parameters where not parametrized on multi-eGO.\n; If this is not expected, check the reference topology."
    else:
        df.rename(columns={df.columns[0]: f"; {df.columns[0]}"}, inplace=True)
        if len(col_space) == 0: return df.to_string(index=False, header=header)
        else:                   return df.to_string(index=False, header=header, col_space = col_space)

def define_c12(sigma, epsilon):
    kT = 2.49 #kJ/mol
    #solving LJ equation = KT and set it equals to a c12 repulsion
    appo = 2*np.power(sigma, 6)/kT * (-epsilon + np.sqrt(epsilon*(epsilon + kT)))
    return np.power(appo, 2)*kT

def c12_pairs_comb_rule(pairs_noH,atoms_types_noH,PAIRS_SCALE_FACTOR):
    
    c12_i = []
    c12_j = []

    for i in range(len(pairs_noH)):
        ai = pairs_noH["ai"].to_numpy()[i]
        aj = pairs_noH["aj"].to_numpy()[i]
        c12_i.append(atoms_types_noH["c12"].loc[atoms_types_noH["atom"]==ai])
        c12_j.append(atoms_types_noH["c12"].loc[atoms_types_noH["atom"]==aj])
    
    c12_i = np.array(c12_i)
    c12_j = np.array(c12_j)

    pairs_noH["c6"]  = 0.0
    pairs_noH["c12"] = np.sqrt(c12_i*c12_j)*PAIRS_SCALE_FACTOR

    return pairs_noH

def get_masses(bonds,atoms_types_noH):
    
    mass_i = []
    mass_j = []
    for i in range(len(bonds)):
        ai = bonds["ai"].to_numpy()[i]
        aj = bonds["aj"].to_numpy()[i]
        mass_i.append(float(atoms_types_noH["mass"].loc[atoms_types_noH["atom"]==ai]))
        mass_j.append(float(atoms_types_noH["mass"].loc[atoms_types_noH["atom"]==aj]))
    
    mass_i = np.array(mass_i)
    mass_j = np.array(mass_j)
    reduced_mass = mass_i*mass_j/(mass_j+mass_i)
    
    return reduced_mass

def save_gro(gro, mapping_dict, h_atom_num, atom_names_mapped, args):

    gro_pd = pd.DataFrame()
    gro_pd["res"]  = [str(atom.residue.idx+1) + str(atom.residue.name)  for atom in gro.atoms]
    gro_pd["atom"] = [str(atom.name) for atom in gro.atoms]
    gro_pd["ai"]   = [atom.idx+1 for atom in gro.atoms]
    gro_pd["x"]    = [coordinate[0]*0.1 for coordinate in gro.coordinates]
    gro_pd["y"]    = [coordinate[1]*0.1 for coordinate in gro.coordinates]
    gro_pd["z"]    = [coordinate[2]*0.1 for coordinate in gro.coordinates]
  
    mask = gro_pd['ai'].isin(h_atom_num)
    gro_noH = gro_pd.copy()
    gro_noH = gro_noH[~mask]
    mapped_indices = [mapping_dict[index] for index in gro_noH["ai"]]
    gro_noH["ai"] = mapped_indices
    gro_noH["atom"] = atom_names_mapped

    with open(f"output/{args.system}/conf_noH.gro", "w") as file:

        file.write(f"{gro.atoms[0].residue.name} no H; {date.strftime('%c')}\n")

        file.write(f"   {len(gro_noH)}\n")

        file.write(f"{dataframe_to_write(gro_noH, False, [8,6,4,7,7,7])}\n")
        file.write(f"   5.00000  5.00000  5.00000")

def save_top(args, atoms_types_noH, atoms_noH, bonds_noH, pairs_noH, angles_noH, dihedrals_noH, impropers_noH, MOL_NAME):
    with open(f"output/{args.system}/topol_noH.top", "w") as file:

        file.write(make_header(args))

        file.write("[ atomtypes ]\n")
        file.write(f"{dataframe_to_write(atoms_types_noH)}\n\n")

        file.write("[ moleculetype ]\n")
        file.write("; Name            nrexcl\n")
        file.write(f"  {MOL_NAME}             3\n\n")

        file.write("[ atoms ]\n")
        file.write(f"{dataframe_to_write(atoms_noH)}\n\n")

        file.write("[ bonds ]\n")
        file.write(f"{dataframe_to_write(bonds_noH)}\n\n")

        file.write("[ pairs ]\n")
        file.write(f"{dataframe_to_write(pairs_noH)}\n\n")

        file.write("[ angles ]\n")
        file.write(f"{dataframe_to_write(angles_noH)}\n\n")

        file.write("[ dihedrals ]\n")
        file.write(f"{dataframe_to_write(dihedrals_noH)}\n\n")

        file.write("[ impropers ]\n")
        file.write(f"{dataframe_to_write(impropers_noH)}\n\n")


def main():

    parser = argparse.ArgumentParser(description="Process input and output files with scaling factors.")
    parser.add_argument("--system", type=str, required=True,  help="input:name of system in input directory")
    parser.add_argument("--bond_scale", type=float, default=1.0, help="Scaling factor for bonds")
    parser.add_argument("--dihedral_scale", type=float, default=1.0, help="Scaling factor for dihedrals")
    parser.add_argument("--angle_scale", type=float, default=1.0, help="Scaling factor for angles")
    parser.add_argument("--impropers_scale", type=float, default=1.0, help="Scaling factor for impropers")
    parser.add_argument("--scale_bonds",  type=bool, default=True, help="Checks frequencies of bond vibration and scale them if necessary")
    parser.add_argument("--scale_masses", type=bool, default=False, help="Rescale masses of heavy atoms after removal of Hs")
    
    args = parser.parse_args()
    if args.scale_masses:
        print("""
              
Scaling masses = True
#WARNING: automatic rescaling works if H atoms are after the connected heavy atom.
If this is not the case modify them by hand
              
              """)
    if args.scale_bonds:
        print("""
Scaling bonds = True
Check frequencies 
              
              """)
    input_files = os.listdir(f"input/{args.system}")
    gro_appo = [f for f in input_files if ".gro" in f]
    top_appo = [f for f in input_files if ".top" in f]

    if not os.path.exists(f"output/{args.system}"):
        os.makedirs(f"output/{args.system}")

    if len(top_appo) != 1:
        print("Input directory contains more or less than 1 .top file! check your input file and select the proper one or add one")
        exit()
    else:
        top = pmd.load_file(f"input/{args.system}/{top_appo[0]}")

    MOL_NAME = top.residues[0].name

    if len(gro_appo) > 1:
        print("Input directory contains more than 1 .gro file! check your input file and select the proper one")
        exit()

    elif len(gro_appo) == 1:
        gro = pmd.load_file(f"input/{args.system}/{gro_appo[0]}")

    #read the topology and extract the needed information
    PAIRS_SCALE_FACTOR = top.defaults.fudgeLJ
    atoms_types = get_atoms_types(top.atoms)
    atoms       = get_atoms(top.atoms)
    bonds       = get_bonds(top.bonds)
    pairs       = get_pairs(top.adjusts)
    angles      = get_angles(top.angles)
    dihedrals   = get_dihedrals(top.dihedrals)
    if len(top.impropers)>0: impropers = get_dihedrals(top.impropers)
    else: impropers = []

    #get hydrogen mapping objects
    mapping_dict, h_atom_num, atoms_types_noH = get_H_mapping(atoms_types, args)
    mapping_dict, h_atom_num, atoms_noH       = get_H_mapping(atoms, args)

    #reindex the new dataframe
    atoms_noH["nr"]         = np.arange(1,len(atoms_noH)+1, 1)
    atoms_noH["atom"]       = np.arange(1,len(atoms_noH)+1, 1)
    atoms_types_noH["atom"] = np.arange(1,len(atoms_types_noH)+1, 1)

    #define c6 and c12 in mego fashion
    atoms_types_noH["c6"] = 0.0
    atoms_types_noH["c12"] = define_c12(atoms_types_noH["sigma"], atoms_types_noH["epsilon"])

    #For each entry of the topology copy the dataframe, remove all entries with H and then
    # map the atom index to the new indexing without H
    bonds_noH     = convert_bonds(bonds, h_atom_num, mapping_dict)
    pairs_noH     = convert_pairs(pairs, h_atom_num, mapping_dict)
    angles_noH    = convert_angles(angles, h_atom_num, mapping_dict)
    dihedrals_noH = convert_dihedrals(dihedrals, h_atom_num, mapping_dict)
    if len(impropers)>0: impropers_noH = convert_impropers(impropers, h_atom_num, mapping_dict)
    else: impropers_noH = pd.DataFrame()

    #calculate oscillation fewquencies of bonds 
    reduced_masses = get_masses(bonds_noH, atoms_types_noH)
    f = np.sqrt(bonds_noH["c1"].to_numpy()/reduced_masses)
    
    #TODO compare this frequencies to the dt step and scale them down if necessary
    GROM_THRESHOLD = 10
    dt_mego = 5e-3

    #print("T = ",np.pi*2/f)
    if np.any(np.pi/(f*dt_mego)<GROM_THRESHOLD) and args.scale_bonds:
        print("""###FAST OSCILLATION FOUND --> RESCALE BONDS""")

        f_max = np.max(f)
        new_bonds = bonds_noH["c1"].to_numpy() / (f_max * dt_mego /(2*np.pi) * GROM_THRESHOLD* 1.01 )**2
        print("Rescaled all bonds")
        ff = np.sqrt(new_bonds/reduced_masses)
    else:
        new_bonds = bonds_noH["c1"].to_numpy() 
        ff = np.sqrt(new_bonds/reduced_masses)

    bonds_noH["c1"] = new_bonds
    #print(bonds_noH)
    #insert pairs c12 based on the combination rule of mego c12 with fudgeLJ
    pairs_noH = c12_pairs_comb_rule(pairs_noH, atoms_types_noH, PAIRS_SCALE_FACTOR)

    #remove unnecessary columns
    pairs_noH = pairs_noH.drop(columns=['type'])
    atoms_types_noH = atoms_types_noH.drop(columns=['atom','sigma', 'epsilon', 'name'])

    #translate names to mego standard
    atom_names_mapped = []
    for name in atoms_noH["name"] :
        if name in list(from_ff_to_multiego):
            atom_names_mapped.append(from_ff_to_multiego[name])
        else: 
            atom_names_mapped.append(name)
    atoms_noH["name"] = atom_names_mapped
    #atoms_noH["type"] = atom_names_mapped
    #atoms_types_noH["name"] = atom_names_mapped
    atoms_types_noH["c12"] = ['{:1e}'.format(item) for item in atoms_types_noH["c12"]]
    atoms_types_noH = atoms_types_noH.drop_duplicates()

    #write gro file with consistent atom names
    if len(gro_appo) > 0: save_gro(gro, mapping_dict, h_atom_num, atom_names_mapped, args)

    # write output topology
    save_top(args, atoms_types_noH, atoms_noH, bonds_noH, pairs_noH, angles_noH, dihedrals_noH, impropers_noH, MOL_NAME)
   
    atoms_types_noH.to_csv(f'output/{args.system}/custom_c12.csv', index=False)  
    print("""
\n\nFinished converting the topology.
Don't forget to fix the dihedrals!! ;-)
          """)

if __name__ == "__main__":
    main()
