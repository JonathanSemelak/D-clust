import argparse
import sys
import math
import numpy as np

# def calc_dihedral(data):
#     ''' The order of the input elements is the natural definition.
#     A --> B --> C --> D
#     '''
#     v1 = data[:,1,:] - data[:,0,:]
#     v2 = data[:,2,:] - data[:,1,:]
#     v3 = data[:,3,:] - data[:,2,:]
#
#     n1 = np.cross(v1,v2)
#     n2 = np.cross(v2,v3)
#
#     dot = (n1 * n2).sum(axis=1)
#     norm1 = np.linalg.norm(n1,axis=1)
#     norm2 = np.linalg.norm(n2,axis=1)
#
#     phi = np.arccos(dot / (norm1 * norm2))
#     return np.degrees(phi)
#
# import numpy as np

def calculate_dihedral(trajectory, i, j, k, l):
    dihedrals = []
    for step in trajectory:
        A = step[i]
        B = step[j]
        C = step[k]
        D = step[l]

        v1 = B - A
        v2 = C - B
        v3 = D - C

        n1 = np.cross(v1, v2)
        n2 = np.cross(v2, v3)

        n1_mag = np.linalg.norm(n1)
        n2_mag = np.linalg.norm(n2)

        cos_phi = np.dot(n1, n2) / (n1_mag * n2_mag)
        phi = np.arccos(np.clip(cos_phi, -1.0, 1.0)) # Clip for numerical stability

        phi_deg = np.degrees(phi)
        dihedrals.append(phi_deg+180.0)

    return np.array(dihedrals)

script_description = """
    ************************************************************
    *                                                          *
    *           Welcome to the D-clust Python Script!          *
    *                                                          *
    *   This tool will help use to extract temporal courses    *
    *   of dihedral angles from your trajectory file and to    *
    *  analize them using scipy, sklearn and DADApy libraries. *
    *                                                          *
    *                                                          *
    *                                                          *
    ************************************************************
"""

def print_welcome_message(script_description):
    print(script_description)

#-------------------------------------------
# Create the parser
parser = argparse.ArgumentParser(description=script_description,formatter_class=argparse.RawDescriptionHelpFormatter)
# Define the command-line arguments
parser.add_argument('-i', '--input', required=True, help='Input file name')
parser.add_argument('-o', '--output', required=True, help='Output file name')
parser.add_argument('-d', '--dihelist', default='none', help="A text file with the atom index of each dihedral to be extracted (not needed if format is 'dihe')")
parser.add_argument('-f', '--format', required=True, help="Input file format ('xyz', 'netcdf'  or 'dihe')")

# If no arguments are provided, print the description and exit
if len(sys.argv) == 1:
    parser.print_help()
    sys.exit()

# Parse the arguments
args = parser.parse_args()

# Assign values from args
input_name = args.input
output_name = args.output
dihelist_name = args.dihelist
file_format = args.format

# Conditional import based on file extension
if (file_format=='xyz'): from ase.io import read
if (file_format=='netcdf'): from scipy.io import netcdf_file

# Call the function at the beginning of your main script execution
if __name__ == "__main__":
    print_welcome_message(script_description)
    # Rest of your script follows here...

# The main code starts here <------------------------------------------------------------------

# Reads data
if (file_format=='xyz'):  # XYZ file case
    print("\nCoordinates from the xyz file will be read using the ASE library\n")
    print("\nReading file...\n")
    trajectory = read(input_name, index=':')
    nsteps = len(trajectory)
    natoms = len(trajectory[0])
    coordinates = np.empty((nsteps, natoms, 3))
    for i, frame in enumerate(trajectory):
        coordinates[i] = frame.get_positions()
elif (file_format=='netcdf'): # NETCDF file case
    print("\nCoordinates from the netcdf file will be read using the scipy library\n")
    print("\nReading file...\n")
    trajectory = netcdf_file(input_name, 'r')
    coordinates = np.array(trajectory.variables['coordinates'].data)
    nsteps = len(coordinates)
    natoms = len(coordinates[0])

if ((file_format=='xyz') or (file_format=='netcdf')): #In this case a file with the dihe definition must be provided
    dihelist=np.loadtxt(dihelist_name,dtype='int')
    if (len(np.shape(dihelist))==1): #1D array recieved
        old_dihelist=dihelist
        dihelist   = np.empty((len(old_dihelist)-3, 4))
        for i in range(0,len(old_dihelist)-3):
            dihelist[i]=[old_dihelist[i],old_dihelist[i+1],old_dihelist[i+2],old_dihelist[i+3]]
        dihelist=dihelist.astype(int)
    ndihe = len(dihelist)
    dihetraj = np.empty((nsteps,ndihe))
    print("\nCalculating dihedrals temporal traces...\n")
    for i in range(0,ndihe):
        print("--> Dihedral", i+1, "(out of ",ndihe, ")")
        dihetraj[:,i]=calculate_dihedral(coordinates,dihelist[i][0],dihelist[i][1],dihelist[i][2],dihelist[i][3])
    print("\nThis results will be saved to 'dihetraj.dat' file...\n")
    fmt = ['%d'] + ['%.4f'] * ndihe
    indexes=np.arange(nsteps)
    np.savetxt('dihetraj.dat',np.column_stack((indexes,dihetraj)),fmt=fmt)
else: # dihe case (a dihetraj file is provided
    # Open the file and read one line to determine the number of columns (dihe+1)
    print("\nDihedrals time evolution will be read directly from input file\n")
    print("\nReading file...\n")
    with open(input_name, 'r') as file:
      first_line = file.readline()
    ndihe = len(first_line.split())-1
    # Load the data, skipping the first column
    dihetraj=np.loadtxt(input_name,usecols=range(1, ndihe+1))
    nsteps=len(dihetraj)

#From now onwards, we will be working with the dihetraj array (nsteps,ndihe)
from dadapy import Data
from dadapy import plot as pl
import matplotlib.pyplot as plt

dihetraj=dihetraj[::10]
dihetraj = np.clip(dihetraj, 0.001, 359.999) # Clip for numerical stability
dihetraj=dihetraj*np.pi/180.0 #Converts to radians

print(np.max(dihetraj/(2*np.pi)))

# initialise a Data object
d_dihedrals = Data(dihetraj, verbose=False)
# compute distances by setting the correct period
d_dihedrals.compute_distances(maxk=dihetraj.shape[0]-1, period=2.*np.pi)
# estimate the intrinsic dimension
d_dihedrals.compute_id_2NN()
# ID scaling analysig using two different methods
ids_2nn, errs_2nn, scales_2nn = d_dihedrals.return_id_scaling_2NN()
ids_gride, errs_gride, scales_gride = d_dihedrals.return_id_scaling_gride(range_max=1024)
col = 'darkorange'
plt.plot(scales_2nn, ids_2nn, alpha=0.85)
plt.errorbar(scales_2nn, ids_2nn, errs_2nn, fmt='None')
plt.scatter(scales_2nn, ids_2nn, edgecolors='k',s=50,label='2nn decimation')
plt.plot(scales_gride, ids_gride, alpha=0.85, color=col)
plt.errorbar(scales_gride, ids_gride, errs_gride, fmt='None',color=col)
plt.scatter(scales_gride, ids_gride, edgecolors='k',color=col,s=50,label='2nn gride')
plt.xlabel(r'Scale',size=15)
plt.ylabel('Estimated ID',size=15)
plt.xticks(size=15)
plt.yticks(size=15)
plt.legend(frameon=False,fontsize=14)
plt.tight_layout()
plt.show()
