
import argparse
import sys
import numpy as np
from parmed.amber import NetCDFTraj

script_description = """
    *************************************************************
    *                                                           *
    *     Welcome to get_frames_from_list Python Script!        *
    *                                                           *
    *   This script takes a list with frame numbers, a          *
    *   trajectory file in netcdf format, and extracts the      *
    *   frames from the list, generating a new trajectory file. *
    *    It requires numpy and ParmED                           *
    *                                                           *
    *************************************************************
"""
#-------------------------------------------
# Create the parser
parser = argparse.ArgumentParser(description=script_description,formatter_class=argparse.RawDescriptionHelpFormatter)
# Define the command-line arguments
parser.add_argument('-i', '--input', required=True, help='Input trajectory file name')
parser.add_argument('-l', '--list', required=True, help='List with frames (single column, integers)')
parser.add_argument('-o', '--output', required=True, help='Output trajectory file name')

# If no arguments are provided, print the description and exit
if len(sys.argv) == 1:
    parser.print_help()
    sys.exit()

# Parse the arguments
args = parser.parse_args()

input_traj_name=args.input
input_list_name=args.list
output_traj_name=args.output

print("\n Coordinates from the netcdf file will be read using the ParmED library\n")
print("\n Reading file "+ input_traj_name +" ...\n")

trajectory = NetCDFTraj.open_old(input_traj_name)
coordinates = np.array(trajectory.coordinates)
natoms = len(coordinates[0])
selected_frames = np.loadtxt(input_list_name,dtype='int')


# Open the NetCDF trajectory file for writing
crd = NetCDFTraj.open_new(output_traj_name, natom=natoms, box=True, crds=True, vels=False, frcs=False)
box=trajectory.box

selected_coordinates=coordinates[selected_frames]

print("\n Processing ...\n")

for j in range(0,len(selected_coordinates)):
    crd.add_coordinates(selected_coordinates[j])
    crd.add_box(box[j])

crd.close()
print("\n --> Frames listed in "+input_list_name+" saved to "+output_traj_name+ " file \n")
