## Convert topology

The script is ment to generate a topology for multi-eGO starting from the one of the target simulation.
The idea is to: 
- remove all hydrogens 
- reindexing all indeces according to the new enumeration
- translate atom names to a multi-eGO consistent version
- (eventually) rescale all bonds, dihedrals, angles forces to be compatible with the 5fs timestep and the new masses

### What to do before using the converter:

Before running the script the user should rename the atom types accordingly to multi-eGO type. In doing so, one should also modify the masses of the atoms in which the hydrogen was removed.

- Make a copy of the starting topology top_copy.top
- modify the the masses in top_copy.top
- run the script using as input the modified top_copy.top

### OUTPUT

- mego consistent topology
- mego gro with correct atom names
- custom dictionary to pass to multiego to consider new c12s

# How to run
python convert_top.py --system system_name

system_name is the directory name in input directory
the output is stored in output/system_name/