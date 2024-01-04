# D-clust (using DADApy)

## Welcome to the D-clust Python Script!  
This code is pretty much an automatization tool that follows the DADApy 
official tutorial, but it makes it easier for AMBER users because AMBER 
netcdf files can be provided as input directly.

It can extract temporal courses of dihedral angles from your
trajectory file and to analize them with the DADApy library using the 
Advanced Density Peaks (ADP) algorithm. It also recieves a file with suche
time courses if you prefer so.

It is recommended to check the DadaPy and ADP papers:

- Glielmo, A., Macocco, I., ... & Laio, A. (2022). DADApy: Distance-based
analysis of data-manifolds in Python. Patterns, 3(10).

- Rodriguez, A., & Laio, A. (2014). Clustering by fast search and find of
density peaks. Science, 344(6191), 1492-1496.

Also, make sure to check the DADApy library repo:

https://github.com/sissa-data-science/DADApy

## Usage
Run the script from the command line with the required arguments.

```console
python D-clust.py [-h] -i INPUT -f FORMAT [-d DIHELIST] [other optional arguments]
```

### Command-Line Arguments
- `-i`, `--input`: Input file name (required).
- `-d`, `--dihelist`: A text file with the atom index of each dihedral to be extracted (required unless the format is 'dihe').
- `-f`, `--format`: Input file format ('xyz', 'netcdf' or 'dihe', required). 
- `-id`, `--id`: Intrinsic dimension (integer, if not provided, it will be estimated).
- `-v`, `--visualize`: Plot on screen the intrinsic dimension scaling ('True' or 'False', if 'False', it only saves it as a .svg file, default: 'False').
- `-ha`,`--halo`: Use halo for ADP ('True' or 'False', default: 'False').
- `-z`, `--zvalue`: Z value for ADP (float, default: 3.5).
- `wt`, `--writetrajs`: Write a trajectory file for each cluster ('True' or 'False', default: 'False').
- `-wf`,`--writefreq`: Writting frequence (for --writetrajs/-wt option, intgeter, default: 1).
- `-nj`,`--njobs`: Number of threads for ADP calculation (integer, default: 1).
- `-s` ,`--slice`: Analize a slice of the data (two integers, not used by default, frame count starts with zero).
- `-rc`, `--randomchoice`: Make a random choice of --randomchoice/-rc frames (integer, not used by default).
- `-h`, `--help`: Show help message and exit.

### Comments about "dihelist" and the "dihe" format

The `dihelist` file can be provided in two formats: A single-column file, containing the number of the atoms forming a "chain" of dihedral angles, or a four-columng file, in which the mask of each dihedral is specified. In every case, the atom count starts with 0.

Single-column example:
```console
5
7
8
9
4
```

Four-column example:
```console
5 7 8 9
7 8 9 4
```
In both cases, the dihedral angles that will be considered are: 5-7-8-9 and 7-8-9-4.

If the `format` is 'dihe', a file with NDIHE+1 columns should be provided, containing tha time evolution of NDIHE dihedral angles (columns 2 to NDIHE+1) and a column containing the frame number.

Example:

```console
1 dihe_1_value_in_frame_1 dihe_1_value_in_frame_1 dihe_3_value_in_frame_1 ...
2 dihe_1_value_in_frame_2 dihe_1_value_in_frame_2 dihe_3_value_in_frame_2 ...
3 dihe_1_value_in_frame_3 dihe_1_value_in_frame_3 dihe_3_value_in_frame_3 ...
.
.
.
```
This option is usefull for saving-time purposes only. Also, a file in this format will be generated as output (with the name 'dihetraj.dat') if `format` is 'netcdf' of `xyz`.


## Acknowledgement
Please cite the corresponding DADApy and ADP papers in case you use this tool.


