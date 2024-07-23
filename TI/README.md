# Description
The script should read a json file or a argparse with the atom indeces corresponding to the six dof for the restraint 

The atom_idx will be read corresponding to the following "structure"

    protein_c,protein_b,protein_a,ligand_A,lig_B,lig_C

where the interactions are

- bond  : aA
- angle1: baA
- angle2: aAB
- dih1  : cbaA
- dih2  : baAb
- dih3  : aABC

The script will print the lines which should be appended in the topology for TI

The x_0 in the json file will be read in this order:

- bond_0
- angle1_0
- angle2_0
- dih1_0
- dih2_0
- dih3_0


### Example of json file
```
{
    "atom_idx": [1, 2, 3, 4, 5, 6],          #c,b,a,A,B,C
    "x_0"     : [0.5, 100, 90, 170, 89, 90]  #d_0, a1_0, a2_0, dih1_0, dih2_0, dih3_0
    "k"       : [4184, 41.84, 41.84, 41.84, 41.84, 41.84]   #associated harmonic constants
}
```

# Calculate barint
```
mkdir all_md ;for j in {0..25}; do echo ${j}; cp Lambda_${j}/MD/md${j}.xvg all_md; done; cd all_md ; gmx_mpi bar -f md* -o -oi; cd ..
```
or 
```
gmx_mpi bar -f md* -o -oi
```

# Generate mdps
Go in mdp directory (e.g. mdps/MD) and run:
```
perl write_mdp.pl file.mdp
```

