# Master thesis

Automatizing MDpocket for the detection of allosteric pockets and its descriptors in GPCR dynamics simulations.

## Instructions

This guide describes how to use the automatization of MDpocket for pocket detection over protein simulations.

### Prerequisites

* [AmberTools](https://ambermd.org/AmberTools.php)
* [GROMACS/5.1.2](https://manual.gromacs.org/documentation/5.1.2/install-guide/)
* [fpocket](http://sourceforge.net/projects/fpocket). This will provide you a file like fpocket-src-1.0.tgz.
In order to install fpocket now, use the following series of commands in a command line:
```console
tar -xzf fpocket-src-2.0.tgz
cd fpocket-src-2.0/
make
make test
```
If the make and make test command yield no errors, your installation can be completed by
typing (if you have administrators rights :
```console
sudo make install
```

### MDPOCKET
MDpocket is an open-source tool for tracking small molecule binding sites on molecular dynamics (MDs) trajectories. 
It is based on fpocket  geometry-based cavity detection algorithm. 
In order to run MDpocket, you should have at minimum a Pentium III 500 Mhz with 128Mb of
RAM.   This   program   was   co­developed   and   tested   under   the   following   Linux   distributions   :
openSuse 10.3 (and newer), Centos 5.2, Fedora Core 7, Ubuntu 8.10 as well as Mac OS X (10.5 &
10.6). Make sure you have gcc installed.

### ALGORITHM INPUT

As algorithm input we need a directory path where all the trajectories files along with it's model file is found.
As it has to be implemented into GPCRmd database, the algorithm automatically matches each trajectory with its model file
thanks to the nomenclature of GPCRmd files.

Example:

We can download the trajectories and the pdb files of the Beta-2 adrenergic receptor in complex with epinephrine from [here](https://submission.gpcrmd.org/dynadb/dynamics/id/117/).
We save the files in a directory. We can see that all those files are have the same format:<br />
[trj/dyn]_[dynamicsID].extension <br />
The algorithm matches the pdb file with its trajectory using the dynamicsID. So to run the algorithm we just need to execute the following command:

```console
python automatization_MDpocket [path_to_directory]
```
where path_to directory is the path where all the trajectories and its model files are stored.

### WORKFLOW

* Matching pdb files with trajectories: First of all as mentioned before the algorithm matches each trajectory with its corresponding model file (pdb).

* Creating MDpocket input: MDpocket uses a list of pdb files as input. For this reason each frame of the trajectory is converted into a pdb file. We will create a folder named snapshots_input_NameTrajectoryFile where all pdb files will be stored
  with the name: md + number of the snapshot in the simulation + .pdb. For a matter of computation time it will just be created a 1 pdb file every 5 frames but it can be changed with the variable pdb_step.
  The pdb files created will only contain the protein. Finally, it will be created a file (snapshot_list_file_TrajectoryFileName) where the path to every pdb file is found. This file will be used as input for MDpocket.
  
* Running MDpocket: MDpocket is run over the simulation and its output is stored in the folder mdpocket. To see more information about the output of MDpocket you can go [here](http://fpocket.sourceforge.net/manual_fpocket2.pdf) (page 21-24).

* Extracting pocket coordinates: from the output file of MDpocket mdpout_freq.dx which has a grid that contains a measure of frequency of how many times the pocket was open during a MD trajectory. The frequency of the pocket is expressed as the number of Voroni Vertices in a 8 Armstrong cube around each grid point per snapshot, so the more a cavity is conserved the highers is the value. In order to detect transient pockets it is extracted the grid points with an isovalue equal or higher than 3. As a result, we obtain all the coordinates detected as a pockets with an isovalue equal or higher than 3. These coordinates can be visualized with the following file: Extract_ISO_output_3_TrajectoryFileName.
  
* Clustering coordinates into pockets: Using density-based spatial clustering of applications with noise (DBSCAN) algorithm the coordinates are clustered into different pockets. The idea behind DBSCAN is that a point belongs to a cluster if it is close to many points from that cluster. In addition, it is albe to detect outliers (coordinates that are far from the others) allowing not to take into consideration coordinate points that are outliers.
