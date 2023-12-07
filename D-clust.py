import argparse
import sys
import math
import numpy as np
import warnings
from plot_functions import *

# Suppress specific UserWarnings from dadapy
warnings.filterwarnings('ignore', message="data type is float64: most methods work only with float-type inputs", category=UserWarning, module='dadapy')

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
    *   This tool will help you to extract temporal courses    *
    *   of dihedral angles from your trajectory file and to    *
    *   analize them with the DADApy library.                  *
    *                                                          *
    *  It's pretty much an automatization tool that follows    *
    *  the DADApy official tutorial, but it makes it easier    *
    *  for AMBER users.                                        *
    *                                                          *
    *                                                          *
    ************************************************************
"""

def print_welcome_message(script_description):
    print(script_description)

def check_bool(inputvariable,rightvariable):
    if inputvariable not in ['True','False']:
        sys.exit("Error: "+ rightvariable + " must be 'True' or 'False'")

#-------------------------------------------
# Create the parser
parser = argparse.ArgumentParser(description=script_description,formatter_class=argparse.RawDescriptionHelpFormatter)
# Define the command-line arguments
parser.add_argument('-i', '--input', required=True, help='Input file name')
parser.add_argument('-d', '--dihelist', default='none', help="A text file with the atom index of each dihedral to be extracted (not needed if format is 'dihe')")
parser.add_argument('-f', '--format', required=True, help="Input file format ('xyz', 'netcdf'  or 'dihe')")
parser.add_argument('-id', '--id', default=0, help="Intrinsic dimension")
parser.add_argument('-v', '--visualize', default="False", help="Intrinsic dimension")
parser.add_argument('-ha', '--halo', default="False", help="Use halo for ADP")
parser.add_argument('-z', '--zvalue',  default=3.5, help="Z value for ADP")

# If no arguments are provided, print the description and exit
if len(sys.argv) == 1:
    parser.print_help()
    sys.exit()

# Parse the arguments
args = parser.parse_args()

# Assign values from args
input_name = args.input
dihelist_name = args.dihelist
file_format = args.format
z_value = args.zvalue
ID = args.id

# Checks the variables are str True or False before converting to bool
check_bool(args.visualize,"--visualize ( -v)")
visualize = args.visualize == "True"

check_bool(args.halo,"--halo ( -ha)")
halo = args.halo == "True"

# Conditional import based on file extension
if (file_format=='xyz'): from ase.io import read
if (file_format=='netcdf'): from parmed.amber import NetCDFTraj



# Call the function at the beginning of your main script execution
if __name__ == "__main__":
    print_welcome_message(script_description)
    # Rest of your script follows here...

# The main code starts here <------------------------------------------------------------------

# Reads data
if (file_format=='xyz'):  # XYZ file case
    print("\n Coordinates from the xyz file will be read using the ASE library\n")
    print("\n Reading file...\n")
    trajectory = read(input_name, index=':')
    nsteps = len(trajectory)
    natoms = len(trajectory[0])
    coordinates = np.empty((nsteps, natoms, 3))
    for i, frame in enumerate(trajectory):
        coordinates[i] = frame.get_positions()
elif (file_format=='netcdf'): # NETCDF file case
    print("\n Coordinates from the netcdf file will be read using the ParmED library\n")
    print("\n Reading file...\n")
    trajectory = NetCDFTraj.open_old(input_name)
    coordinates = np.array(trajectory.coordinates)
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
    print("\n Calculating dihedrals temporal traces...\n")
    for i in range(0,ndihe):
        print(" --> Dihedral", i+1, "( out of ",ndihe,")")
        dihetraj[:,i]=calculate_dihedral(coordinates,dihelist[i][0],dihelist[i][1],dihelist[i][2],dihelist[i][3])
    print("\n These results will be saved to 'dihetraj.dat' file...\n")
    fmt = ['%d'] + ['%.4f'] * ndihe
    indices=np.arange(nsteps)
    np.savetxt('dihetraj.dat',np.column_stack((indices,dihetraj)),fmt=fmt)
else: # dihe case (a dihetraj file is provided
    # Open the file and read one line to determine the number of columns (dihe+1)
    print("\n Dihedrals time evolution will be read directly from input file\n")
    print("\n Reading file...\n")
    with open(input_name, 'r') as file:
      first_line = file.readline()
    ndihe = len(first_line.split())-1
    # Load the data, skipping the first column
    dihetraj=np.loadtxt(input_name,usecols=range(1, ndihe+1))
    nsteps=len(dihetraj)

#From now onwards, we will be working with the dihetraj array (nsteps,ndihe)
from dadapy import Data
# from dadapy import plot as pl
import matplotlib.pyplot as plt

dihetraj = np.clip(dihetraj, 0.001, 359.999) # Clip for numerical stability
dihetraj = dihetraj*np.pi/180.0 #Converts to radians

# initialise a Data object
d_dihedrals = Data(dihetraj, verbose=False)
# compute distances by setting the correct period
d_dihedrals.compute_distances(maxk=dihetraj.shape[0]-1, period=2.*np.pi)

# estimate the intrinsic dimension
if (ID == 0):
    print("\n The scaling of the Intrinsic Dimension will be evaluated using the 2nn and GRIDE methods\n")
    print("\n Computing ID...\n")

    # ID scaling analysig using two different methods
    ids_2nn, errs_2nn, scales_2nn = d_dihedrals.return_id_scaling_2NN()
    ids_gride, errs_gride, scales_gride = d_dihedrals.return_id_scaling_gride(range_max=1024)

    print("\n 2nn ID scaling:\n")
    print("\n Scale  | Estimated ID  | Error on ID:\n")

    for i in range(0, len(ids_2nn)):
        print(f" {scales_2nn[i]:.3f} {ids_2nn[i]:.3f} {errs_2nn[i]:.3f}")

    print("\n GRIDE ID scaling:\n")
    print("\n Scale  | Estimated ID  | Error on ID:\n")

    for i in range(0, len(ids_gride)):
        print(f" {scales_gride[i]:.3f} {ids_gride[i]:.3f} {errs_gride[i]:.3f}")

    if (visualize): #This was taken from the DADApy tutorial
        print("\n Showing plot ID scaling plot...\n")
        plot_ID_scaling(scales_2nn,ids_2nn,errs_2nn,scales_gride,ids_gride, errs_gride, savefig="ID_scaling.svg", showplot=True)
    else:
        plot_ID_scaling(scales_2nn,ids_2nn,errs_2nn,scales_gride,ids_gride, errs_gride, savefig="ID_scaling.svg", showplot=False)
    print("\n ---> ID scaling plot saved to ID_scaling.svg file\n")

    print("\n Assuming a plateu is reached\n")
    print("\n Make sure this is the case by visualizing the ID scaling!\n")
    print("\n The ID will be approximated as the maximum between the minimum ID estimation of 2nn and GRIDE\n")

    ID=int(max(np.min(ids_2nn),np.min(ids_gride)))

    print("\n Estimated ID:", int(ID) )
else:
    ID=int(ID)
    print("\n The Intrinsic Dimension (ID) was given as input:\n")
    print("\n Input ID:", int(ID),"\n" )

print("\n Performing Advanced Density Peaks (ADP) analysis:\n")
print("\n Clusterizing...\n")

# cluster data via Advanced Density Peak
d_dihedrals.set_id(ID)
d_dihedrals.compute_clustering_ADP(Z=float(z_value),halo=halo);
n_clusters = len(d_dihedrals.cluster_centers)

print("\n Performing Advanced Density Peaks (ADP) analysis:\n")

if (halo):
    print("\n Number of clusters found:", int(n_clusters) ,"(Z value =", z_value, " with halo points) \n")
else:
    print("\n Number of clusters found:", int(n_clusters) ,"(Z value =", z_value, " no halo points) \n")

if (visualize):
    print("\n Showing dendrogram plot...\n")
    get_dendrogram_custom(d_dihedrals, cmap='Set2', savefig="dendogram.svg", logscale=False,showplot=True)
else:
    get_dendrogram_custom(d_dihedrals, cmap='Set2', savefig="dendogram.svg", logscale=False,showplot=False)
print("\n ---> Dendrogram plot saved to dendogram.svg file\n")



# Cluster populations
cluster_indices = d_dihedrals.cluster_indices
populations = [ len(el) for r_,el in enumerate(cluster_indices)]

print("\n Clusters population:\n")
print("\n #Cluster  |   #Frames:\n")
for i in range(0, n_clusters):
    print(f" {i:.0f} {populations[i]:.0f}")

# Saves list of indices for each cluster
print("\n Saving clurters indices...\n")
print("\n #Cluster  |   #Frames:\n")
for i in range(0,n_clusters):
    ith_cluster_indices=cluster_indices[i]
    ith_cluster_indices=np.array(ith_cluster_indices)
    ith_cluster_indices_filename='cluster_'+str(int(i))+'_indices.dat'
    np.savetxt(ith_cluster_indices_filename,ith_cluster_indices,fmt='%i')
    print(" --> Indices from cluster #"+str(int(i))+" saved in file "+ith_cluster_indices_filename)


# Cluster centers. In the original trajecotory these frames are given by (center + 400) * 10
print("\n Clusters centers:\n")
print("\n #Cluster  |   Center:\n")

centers=d_dihedrals.cluster_centers
for i in range(0, n_clusters):
    print(f" {i:.0f} {centers[i]:.0f}")

# Write stuff
write_trajs=True
freq_write=1

# Cluster centers. In the original trajecotory these frames are given by (center + 400) * 10
print("\n Saving trajectories for each cluster (writing frequence = ", int(freq_write) ,"):\n")
print("\n Generating trajectory files...\n")

if (write_trajs and not file_format == 'dihe'):
    if (file_format=='netcdf'):
        for i in range(0, n_clusters):
            ith_cluster_indices=cluster_indices[i]
            ith_cluster_indices=np.array(ith_cluster_indices)
            ith_cluster_traj_filename='cluster_'+str(int(i))+'_traj.nc'
            # Open the NetCDF trajectory file for writing
            crd = NetCDFTraj.open_new(ith_cluster_traj_filename, natom=natoms, box=False,
                                     crds=True, vels=False, frcs=False)
            selected_coordinates=coordinates[ith_cluster_indices]
            for j in range(0,len(selected_coordinates)):
                if (j%freq_write == 0):
                    crd.add_coordinates(selected_coordinates[j])
            # Close the file
            crd.close()
            print(" --> Frames belonging to cluster #"+str(int(i))+" saved in trajectory file "+ith_cluster_traj_filename)
        for i in range(0, n_clusters):
            ith_cluster_centroid_filename='cluster_'+str(int(i))+'_center.nc'
            # Open the NetCDF trajectory file for writing
            crd = NetCDFTraj.open_new(ith_cluster_centroid_filename, natom=natoms, box=False,
                                     crds=True, vels=False, frcs=False)
            crd.add_coordinates(coordinates[int(centers[i])])
            crd.close()
            print(" --> Coordinates of the center of cluster #"+str(int(i))+" saved in file "+ith_cluster_centroid_filename)
