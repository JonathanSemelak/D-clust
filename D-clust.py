import argparse
import sys
import math
import numpy as np
import warnings
import matplotlib.pyplot as plt
from matplotlib import cm

def get_dendrogram_custom(Data, cmap="viridis", savefig="", logscale=True, showplot=False):

    """
    This is a modification of the "get_dendrogram" function from DADApy.

    All rights reverved to the DADApy authors.
    .
    .
    .
    Generate a visualisation of the topography computed with ADP.
    The only difference is that it makes optional to show the plot.

    This visualisation fundamentally corresponds to a hierarchy of the clusters built
    with Single Linkage taking as similarity measure the density at the
    border between clusters.
    At difference from classical dendrograms, where all the branches have the same height,
    in this case the height of the branches is proportional to the density of the cluster
    centre.
    To convey more information, the distance in the x-axis between
    clusters is proportional to the population (or its logarithm).

    Args:
        Data: A dadapy data object for which ADP has been already run.
        cmap: (optional) The color map for representing the different clusters,
            the default is "viridis".
        savefig: (str, optional) A string with the name of the file in which the dendrogram
            will be saved. The default is empty, so no file is generated.
        logscale: (bool, optional) Makes the distances in the x-axis between clusters proportional
            to the logarithm of the population of the clusters instead of
            proportional to the population itself. In very unbalanced clusterings,
            it makes the dendrogram more human readable. The default is True.

    Returns:

    """
    # Prepare some auxiliary lists
    e1 = []
    e2 = []
    d12 = []
    L = []
    Li1 = []
    Li2 = []
    Ldis = []
    Fmax = max(Data.log_den)
    Rho_bord_m = np.copy(Data.log_den_bord)
    # Obtain populations of the clusters for fine tunning the x-axis
    pop = np.zeros((Data.N_clusters), dtype=int)
    for i in range(Data.N_clusters):
        pop[i] = len(Data.cluster_indices[i])
        if logscale:
            pop[i] = np.log(pop[i])
    xr = np.sum(pop)
    # Obtain distances in list format from topography
    for i in range(Data.N_clusters - 1):
        for j in range(i + 1, Data.N_clusters):
            dis12 = Fmax - Rho_bord_m[i][j]
            e1.append(i)
            e2.append(j)
            d12.append(dis12)

    # Obtain the dendrogram in form of links
    nlinks = 0
    clnew = Data.N_clusters
    for j in range(Data.N_clusters - 1):
        aa = np.argmin(d12)
        nlinks = nlinks + 1
        L.append(clnew + nlinks)
        Li1.append(e1[aa])
        Li2.append(e2[aa])
        Ldis.append(d12[aa])
        # update distance matrix
        t = 0
        fe = Li1[nlinks - 1]
        fs = Li2[nlinks - 1]
        newname = L[nlinks - 1]
        # list of untouched clusters
        unt = []
        for _ in d12:
            if (e1[t] != fe) & (e1[t] != fs):
                unt.append(e1[t])
            if (e2[t] != fe) & (e2[t] != fs):
                unt.append(e2[t])
            t = t + 1
        myset = set(unt)
        unt = list(myset)
        # Build a new distance matrix
        e1new = []
        e2new = []
        d12new = []
        for j in unt:
            t = 0
            dmin = 9.9e99
            for _ in d12:
                if (e1[t] == j) | (e2[t] == j):
                    if (e1[t] == fe) | (e2[t] == fe) | (e1[t] == fs) | (e2[t] == fs):
                        if d12[t] < dmin:
                            dmin = d12[t]
                t = t + 1
            e1new.append(j)
            e2new.append(newname)
            d12new.append(dmin)

        t = 0
        for _ in d12:
            if (unt.count(e1[t])) & (unt.count(e2[t])):
                e1new.append(e1[t])
                e2new.append(e2[t])
                d12new.append(d12[t])
            t = t + 1

        e1 = e1new
        e2 = e2new
        d12 = d12new

    # Get the order in which the elements should be displayed
    sorted_elements = []
    sorted_elements.append(L[nlinks - 1])

    for jj in range(len(L)):
        j = len(L) - jj - 1
        for i in range(len(sorted_elements)):
            if sorted_elements[i] == L[j]:
                sorted_elements[i] = Li2[j]
                sorted_elements.insert(i, Li1[j])

    add = 0.0
    x = []
    y = []
    label = []
    join_distance = []
    for i in range(len(sorted_elements)):
        label.append(sorted_elements[i])
        j = Data.cluster_centers[label[i]]
        y.append(Data.log_den[j])
        x.append(add + 0.5 * pop[sorted_elements[i]])
        add = add + pop[sorted_elements[i]]
        join_distance.append(add)

    xs = x.copy()
    ys = y.copy()
    labels = label.copy()
    zorder = 0
    for jj in range(len(L)):
        c1 = label.index(Li1[jj])
        c2 = label.index(Li2[jj])
        label.append(L[jj])
        if c1 < len(sorted_elements):
            x.append(join_distance[c1])
        else:
            x.append((x[c1] + x[c2]) / 2.0)
        ynew = Fmax - Ldis[jj]
        y.append(ynew)
        x1 = x[c1]
        y1 = y[c1]
        x2 = x[c2]
        y2 = y[c2]
        zorder = zorder + 1
        plt.plot(
            [x1, x1], [y1, ynew], color="k", linestyle="-", linewidth=2, zorder=zorder
        )
        zorder = zorder + 1
        plt.plot(
            [x2, x2], [y2, ynew], color="k", linestyle="-", linewidth=2, zorder=zorder
        )
        zorder = zorder + 1
        plt.plot(
            [x1, x2], [ynew, ynew], color="k", linestyle="-", linewidth=2, zorder=zorder
        )

    zorder = zorder + 1
    cmal = cm.get_cmap(cmap, Data.N_clusters)
    colors = cmal(np.arange(0, cmal.N))
    plt.scatter(xs, ys, c=labels, s=100, zorder=zorder, cmap=cmap)
    for i in range(Data.N_clusters):
        zorder = zorder + 1
        cc = "k"
        r = colors[labels[i]][0]
        g = colors[labels[i]][1]
        b = colors[labels[i]][2]
        luma = (0.2126 * r + 0.7152 * g + 0.0722 * b) * 255
        if luma < 156:
            cc = "w"
        plt.annotate(
            labels[i],
            (xs[i], ys[i]),
            horizontalalignment="center",
            verticalalignment="center",
            zorder=zorder,
            c=cc,
            weight="bold",
        )
    plt.xlim([-0.02 * xr, xr])
    xname = "Population"
    if logscale:
        xname = "ln(Population)"
        plt.xlim([0, xr])
    plt.xlabel(xname)
    plt.ylabel(r"ln($\rho$)")
    if savefig != "":
        plt.savefig(savefig)
        plt.clf()
    if (showplot):
        plt.show()

def plot_ID_scaling(scales_2nn,ids_2nn,errs_2nn,scales_gride,ids_gride, errs_gride, savefig="", showplot=False):
    col = 'darkorange'
    plt.plot(scales_2nn, ids_2nn, alpha=0.85)
    plt.errorbar(scales_2nn, ids_2nn, errs_2nn, fmt='None')
    plt.scatter(scales_2nn, ids_2nn, edgecolors='k',s=50,label='2nn decimation')
    plt.plot(scales_gride, ids_gride, alpha=0.85, color=col)
    plt.errorbar(scales_gride, ids_gride, errs_gride, fmt='None',color=col)
    plt.scatter(scales_gride, ids_gride, edgecolors='k',color=col,s=50,label='GRIDE')
    plt.xlabel(r'Scale',size=15)
    plt.ylabel('Estimated ID',size=15)
    plt.xticks(size=15)
    plt.yticks(size=15)
    plt.legend(frameon=False,fontsize=14)
    plt.tight_layout()
    if savefig != "":
        plt.savefig(savefig)
        plt.clf()
    if (showplot):
        plt.show()

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
parser.add_argument('-wt', '--writetrajs', default="False", help="Write a trajectory file for each cluster")
parser.add_argument('-wf', '--writefreq', default=1, help="Writting frequence (for --writetrajs/-wt option)")
parser.add_argument('-nj', '--njobs', default=1, help="Number of threads for ADP calculation")
parser.add_argument('-s', '--slice', nargs=2, type=int, default=[0, 0], help='Analize a slice of the data (frame count starts with zero)')
parser.add_argument('-rc', '--randomchoice', default=0, type=int, help="Makes a random choice of --randomchoice/-rc frames")

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
freq_write = int(args.writefreq)
njobs = int(args.njobs)
slice_indices = args.slice
randomchoice = args.randomchoice
use_slice = not (slice_indices[0] == 0 and slice_indices[1] == 0)
take_random = not (randomchoice == 0)

# Checks the variables are str True or False before converting to bool
check_bool(args.visualize,"--visualize ( -v)")
visualize = args.visualize == "True"

check_bool(args.halo,"--halo ( -ha)")
halo = args.halo == "True"

check_bool(args.writetrajs,"--writetrajs ( -wt)")
write_trajs = args.writetrajs == "True"


# Conditional import based on file extension
if (file_format=='xyz'): from ase.io import read
if (file_format=='netcdf'): from parmed.amber import NetCDFTraj

warnings.filterwarnings('ignore', message='Could not find netCDF4 module. Falling back on ', category=UserWarning, module='parmed.amber')


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
    coordinates = np.empty((len(trajectory), len(trajectory[0]), 3))
    for i, frame in enumerate(trajectory):
        coordinates[i] = frame.get_positions()
    if (use_slice):
        print("\n Only a slice will be considered ( from ",slice_indices[0]," to ",slice_indices[1],")\n")
        coordinates=coordinates[slice_indices[0]:slice_indices[1]]
    if (take_random):
        print("\n A random choice of", randomchoice, "frames will be considered\n")
        if (randomchoice > len(coordinates)):
            sys.exit("Error: ", randomchoice ," must be lower than the number of frames to analize")
        indices_array=np.arange(0,len(coordinates))
        indices_randomchoice=np.random.choice(indices_array, randomchoice, replace=False)
        coordinates = coordinates[indices_randomchoice]
        print(" --> The random choice will be saved to random_choice_indices.dat")
        np.savetxt("random_choice_indices.dat",indices_randomchoice,fmt='%i')
    nsteps = len(trajectory)
    natoms = len(trajectory[0])
elif (file_format=='netcdf'): # NETCDF file case
    print("\n Coordinates from the netcdf file will be read using the ParmED library\n")
    print("\n Reading file...\n")
    trajectory = NetCDFTraj.open_old(input_name)
    coordinates = np.array(trajectory.coordinates)
    if (use_slice):
        print("\n Only a slice will be considered ( from ",slice_indices[0]," to ",slice_indices[1],")\n")
        coordinates=coordinates[slice_indices[0]:slice_indices[1]]
    if (take_random):
        print("\n A random choice of", randomchoice, "frames will be considered\n")
        if (randomchoice > len(coordinates)):
            sys.exit("Error: ", randomchoice ," must be lower than the number of frames to analize")
        indices_array=np.arange(0,len(coordinates))
        indices_randomchoice=np.random.choice(indices_array, randomchoice, replace=False)
        coordinates = coordinates[indices_randomchoice]
        print(" --> The random choice will be saved to random_choice_indices.dat")
        np.savetxt("random_choice_indices.dat",indices_randomchoice,fmt='%i')
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
else: # dihe case (a dihetraj file is provided)
    # Open the file and read one line to determine the number of columns (dihe+1)
    print("\n Dihedrals time evolution will be read directly from input file\n")
    print("\n Reading file...\n")
    with open(input_name, 'r') as file:
      first_line = file.readline()
    ndihe = len(first_line.split())-1
    # Load the data, skipping the first column
    dihetraj=np.loadtxt(input_name,usecols=range(1, ndihe+1))
    if (use_slice):
        print("\n Using frames ",slice_indices[0]," to ",slice_indices[1],"\n")
        dihetraj=dihetraj[slice_indices[0]:slice_indices[1]]
    if (take_random):
        print("\n A random choice of", randomchoice, "frames will be considered\n")
        if (randomchoice > len(coordinates)):
            sys.exit("Error: ", randomchoice ," must be lower than the number of frames to analize")
        indices_array=np.arange(0,len(dihetraj))
        indices_randomchoice=np.random.choice(indices_array, randomchoice, replace=False)
        dihetraj = dihetraj[indices_randomchoice]
        print(" --> The random choice will be saved to random_choice_indices.dat")
        np.savetxt("random_choice_indices.dat",indices_randomchoice,fmt='%i')
    nsteps=len(dihetraj)

#From now onwards, we will be working with the dihetraj array (nsteps,ndihe)
from dadapy import Data
# from dadapy import plot as pl
import matplotlib.pyplot as plt

dihetraj = np.clip(dihetraj, 0.001, 359.999) # Clip for numerical stability
dihetraj = dihetraj*np.pi/180.0 #Converts to radians

# initialise a Data object
d_dihedrals = Data(dihetraj, verbose=False,njobs=njobs)
# compute distances by setting the correct period
print("\n Computing distance matrix...\n")
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

print("\n Number of clusters found: ",int(n_clusters)," \n")

print("\n Performing Advanced Density Peaks (ADP) analysis:\n")

if (halo):
    print("\n Number of clusters found:", int(n_clusters) ,"(Z value =", z_value, " with halo points) \n")
else:
    print("\n Number of clusters found:", int(n_clusters) ,"(Z value =", z_value, " without halo points) \n")

if(n_clusters > 1):
    if (visualize):
        print("\n Showing dendrogram plot...\n")
        get_dendrogram_custom(d_dihedrals, cmap='Set2', savefig="dendogram.svg", logscale=False,showplot=True)
    else:
        get_dendrogram_custom(d_dihedrals, cmap='Set2', savefig="dendogram.svg", logscale=False,showplot=False)
    print("\n ---> Dendrogram plot saved to dendogram.svg file\n")
else:
    print("\n Only 1 cluster found, no dendrogram will be generated...\n")

# Prints cluster centers
print("\n Clusters centers:\n")
print("\n #Cluster  |   Center:\n")

centers=d_dihedrals.cluster_centers
for i in range(0, n_clusters):
    print(f" {i:.0f} {centers[i]:.0f}")

# Cluster populations and densities
cluster_indices = d_dihedrals.cluster_indices
populations = [ len(el) for r_,el in enumerate(cluster_indices)]
densities = d_dihedrals.log_den

print("\n Clusters population:\n")
print("\n #Cluster  |   #Frames  | Log Density (center):\n")
for i in range(0, n_clusters):
    print(f" {i:.0f} {populations[i]:.0f} {densities[centers[i]]:.3f}")


# Halo points analysis
if (halo):
    cluster_indices_no_halos=[]
    cluster_indices_halos=[]
    assignment = d_dihedrals.cluster_assignment
    print("\n Performing analysis of halo points...\n")
    for i in range(0,n_clusters):
        ith_cluster_indices=cluster_indices[i]
        ith_cluster_indices_no_halos=[]
        ith_cluster_indices_halos=[]
        for index in ith_cluster_indices:
            ishalo = (assignment[index] == -1)
            if ishalo:
                ith_cluster_indices_halos.append(index)
            else:
                ith_cluster_indices_no_halos.append(index)
        cluster_indices_halos.append(ith_cluster_indices_halos)
        cluster_indices_no_halos.append(ith_cluster_indices_no_halos)

    # Prints halo points
    print("\n Clusters halo points:\n")
    print("\n #Cluster  |   #Halo points\n")
    for i in range(0, n_clusters):
        print(f" {i:.0f} {len(cluster_indices_halos[i]):.0f}")

    # Prints no halo points
    print("\n Clusters no halo points:\n")
    print("\n #Cluster  |   #No halo points\n")
    for i in range(0, n_clusters):
        print(f" {i:.0f} {len(cluster_indices_no_halos[i]):.0f}")


# Saves list of indices for each cluster
print("\n Saving clurters indices...\n")
print("\n #Cluster  |   #Frames:\n")
for i in range(0,n_clusters):
    ith_cluster_indices=cluster_indices[i]
    ith_cluster_indices=np.array(ith_cluster_indices)
    ith_cluster_indices_filename='cluster_'+str(int(i))+'_indices.dat'
    np.savetxt(ith_cluster_indices_filename,ith_cluster_indices,fmt='%i')
    print(" --> Indices from cluster #"+str(int(i))+" saved in file "+ith_cluster_indices_filename)


if (halo):
    for i in range(0,n_clusters):
        ith_cluster_indices_no_halos=cluster_indices_no_halos[i]
        ith_cluster_indices_no_halos=np.array(ith_cluster_indices_no_halos)
        ith_cluster_indices_filename_no_halos='cluster_'+str(int(i))+'_indices_no_halos.dat'
        np.savetxt(ith_cluster_indices_filename_no_halos,ith_cluster_indices_no_halos,fmt='%i')
        print(" --> Indices from cluster #"+str(int(i))+" (no halo points) saved in file "+ith_cluster_indices_filename_no_halos)
    for i in range(0,n_clusters):
        ith_cluster_indices_halos=cluster_indices_halos[i]
        ith_cluster_indices_halos=np.array(ith_cluster_indices_halos)
        ith_cluster_indices_filename_halos='cluster_'+str(int(i))+'_indices_halos.dat'
        np.savetxt(ith_cluster_indices_filename_halos,ith_cluster_indices_halos,fmt='%i')
        print(" --> Indices from cluster #"+str(int(i))+" (halo points) saved in file "+ith_cluster_indices_filename_halos)

# Saves list of indices for each cluster
print("\n Saving clusters assignment...\n")
cluster_assignment_file=np.empty(nsteps)
for i in range(0,nsteps):
     for j in range(0,n_clusters):
          if (i in cluster_indices[j]):
              cluster_assignment_file[i]=j
np.savetxt("cluster_assignment.dat",cluster_assignment_file,fmt='%i')
print(" --> Cluster assignment saved in file cluster_assignment.dat")

# Write trajs
if (write_trajs and not file_format == 'dihe'):
    print("\n Generating trajectory files...\n")
    print("\n Saving trajectories for each cluster (writing frequence = ", int(freq_write) ,"):\n")
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
        if(halo):
            print("\n Saving halo-separeted trajectories for each cluster (writing frequence = ", int(freq_write) ,"):\n")
            for i in range(0, n_clusters):
                ith_cluster_indices_halos=cluster_indices_halos[i]
                ith_cluster_indices_halos=np.array(ith_cluster_indices_halos)
                ith_cluster_traj_filename_halos='cluster_'+str(int(i))+'_halos_traj.nc'
                # Open the NetCDF trajectory file for writing
                crd = NetCDFTraj.open_new(ith_cluster_traj_filename_halos, natom=natoms, box=False,
                                         crds=True, vels=False, frcs=False)
                selected_coordinates=coordinates[ith_cluster_indices_halos]
                for j in range(0,len(selected_coordinates)):
                    if (j%freq_write == 0):
                        crd.add_coordinates(selected_coordinates[j])
                # Close the file
                crd.close()
                print(" --> Frames belonging to cluster #"+str(int(i))+" halo points saved in trajectory file "+ith_cluster_traj_filename_halos)
            for i in range(0, n_clusters):
                ith_cluster_indices_no_halos=cluster_indices_no_halos[i]
                ith_cluster_indices_no_halos=np.array(ith_cluster_indices_no_halos)
                ith_cluster_traj_filename_no_halos='cluster_'+str(int(i))+'_no_halos_traj.nc'
                # Open the NetCDF trajectory file for writing
                crd = NetCDFTraj.open_new(ith_cluster_traj_filename_no_halos, natom=natoms, box=False,
                                         crds=True, vels=False, frcs=False)
                selected_coordinates=coordinates[ith_cluster_indices_no_halos]
                for j in range(0,len(selected_coordinates)):
                    if (j%freq_write == 0):
                        crd.add_coordinates(selected_coordinates[j])
                # Close the file
                crd.close()
                print(" --> Frames belonging to cluster #"+str(int(i))+" without halo points saved in trajectory file "+ith_cluster_traj_filename_no_halos)


print("\n Saving coordinats for each cluster center:\n")
if (file_format=='netcdf'):
    for i in range(0, n_clusters):
        ith_cluster_centroid_filename='cluster_'+str(int(i))+'_center.nc'
        # Open the NetCDF trajectory file for writing
        crd = NetCDFTraj.open_new(ith_cluster_centroid_filename, natom=natoms, box=False,
                                 crds=True, vels=False, frcs=False)
        crd.add_coordinates(coordinates[int(centers[i])])
        crd.close()
        print(" --> Coordinates of the center of cluster #"+str(int(i))+" saved in file "+ith_cluster_centroid_filename)
