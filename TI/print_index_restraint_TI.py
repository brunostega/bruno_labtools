import argparse
import json

def read_json_file(filepath):
    with open(filepath, 'r') as file:
        data = json.load(file)
    atom_idx = data['atom_idx']
    x_0 = data['x_0']
    k = data['k']
    return atom_idx, x_0, k

def main():
    # Set up the argument parser
    parser = argparse.ArgumentParser(description='Generate formatted outputs for molecular interaction setups.')
    parser.add_argument('--atom_idxs', metavar='N', type=int, nargs=6,
                        help='Six integers needed for the formatting.')
    parser.add_argument('--json', type=str, help='json file with data')
    parser.add_argument('--k_bond'    , type=float, default=4184 , help='k bond value')
    parser.add_argument('--k_angle'   , type=float, default=41.84, help='k angle value')
    parser.add_argument('--k_dihedral', type=float, default=41.84, help='k dihedrak value')

    # Parse arguments from command line
    args = parser.parse_args()
    if args.json != None and args.atom_idxs != None:
        print("Either pass the file or the indeces by argparse. Cannot do both")

    if args.atom_idxs != None:# The six integers entered as arguments will be used
        input_numbers = args.atom_idxs

    if  args.json != None:
        input_numbers, x0, k = read_json_file(args.json)


    # Function to format the first output based on specific groupings
    def format_first_output(numbers):
        # Maps of number groups to specific outputs
        groups = [
            [numbers[2], numbers[3]],
            [numbers[1], numbers[2], numbers[3]],
            [numbers[2], numbers[3], numbers[4]],
            [numbers[0], numbers[1], numbers[2], numbers[3]],
            [numbers[2], numbers[3], numbers[4], numbers[5]],
            [numbers[1], numbers[2], numbers[3], numbers[4]]
        ]
        
        output = []
        for group in groups:
            header = "[ " + "_".join(f"a_{num}" for num in group) + " ]"
            numbers_part = " ".join(str(num) for num in group)
            output.append(f"{header}\n {numbers_part}")
        return "\n".join(output)

    # Function to format the second output (hardcoded based on specific numbers)
    def format_second_output(numbers):
        return f"""[ intermolecular_interactions]
    [ bonds ]
    ; ai      aj    type   bA      kA     bB      kB    
    {numbers[2]}   {numbers[3]}     6     {x0[0]}    0.0    {x0[0]}    {k[0]:.2f}

    [ angles ]
    ; ai     aj    ak     type    thA      fcA        thB      fcB
    {numbers[1]}    {numbers[2]}    {numbers[3]}   1       {x0[1]}     0.0       {x0[1]}    {k[1]:.2f}
    {numbers[2]}    {numbers[3]}    {numbers[4]}   1       {x0[2]}     0.0       {x0[2]}    {k[2]:.2f}

    [ dihedrals ]
    ; ai     aj    ak    al    type     thA      fcA       thB      fcB
    {numbers[0]} {numbers[1]} {numbers[2]} {numbers[3]}    2       {x0[3]}    0.0       {x0[3]}     {k[3]:.2f}
    {numbers[1]} {numbers[2]} {numbers[3]} {numbers[4]}    2       {x0[4]}    0.0       {x0[4]}     {k[4]:.2f}
    {numbers[2]} {numbers[3]} {numbers[4]} {numbers[5]}    2       {x0[5]}    0.0       {x0[5]}     {k[5]:.2f}"""

    # Generate and print both outputs
    print(format_first_output(input_numbers))
    print("\n")
    print(format_second_output(input_numbers))

if __name__ == "__main__":
    main()