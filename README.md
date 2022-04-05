# Master thesis

Automatizing MDpocket for the detection of allosteric pockets and its descriptors in GPCR dynamics simulations.

## Instructions

This guide describes how to use the automatization of MDpocket for pocket detection over protein simulations.

### Prerequisites

* [AmberTools](https://ambermd.org/AmberTools.php)
* [GROMACS/5.1.2](https://manual.gromacs.org/documentation/5.1.2/install-guide/)
* [fpocket](http://sourceforge.net/projects/fpocket). This will provide you a file like fpocket-src-1.0.tgz.
* [Miniconda3/4.7.10](https://repo.anaconda.com/miniconda/) 
In order to install fpocket now, use the following series of commands in a command line:
```console
tar -xzf fpocket-src-2.0.tgz
cd fpocket-src-2.0/
make
make test
```
If the make and make test command yield no errors, your installation can be completed by
typing (if you have administrators rights) :
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

### ARGUMENT OPTIONS

The algorithm can work in 2 different ways:

* Automatically detecting the pdb file and the model file from the directory. It can be used if the trajectory file and the pdb file are named using the GPCRmd database format: [trj/dyn]_[dynamicsID].extension. If you want to use the algorithm this way you will need to use the argument -d.
* Parsing the pdb and the trajectory file. If you want to use the algorithm this way you will need to use the arguments -t [trajectory_filename] -p [pdb_filename]

#### REQUIRED ARGUMENTS

##### -d, --input_directory: 

As algorithm input we need a directory path where all the trajectories files along with it's model file is found.
As it has to be implemented into GPCRmd database, the algorithm automatically matches each trajectory with its model file
thanks to the nomenclature of GPCRmd files.

Example:

We can download the trajectories and the pdb files of the Beta-2 adrenergic receptor in complex with epinephrine from [here](https://submission.gpcrmd.org/dynadb/dynamics/id/117/).
We save the files in a directory. We can see that all those files are have the same format:<br />
[trj/dyn]_[dynamicsID].extension <br />
The algorithm matches the pdb file with its trajectory using the dynamicsID. So, to run the algorithm we just need to execute the following command supposing that the pdb files and trajectory files are on the working directory  (the algorithm only works if we execute the algorithm on the current directory):

```console
python automatization_MDpocket -d .
```
where path_to_directory is the path where all the trajectories and its model files are stored. Then the script will automaticaly detect the pockets for the different trajectoris of the Beta-2 adrenergic receptor.

##### -t, --trajectory_file / -p, --pdb_file:
The second manner to run the algorithm is indicating where the trajectory file and the pdb file can be found. The trajectory file should have the xtc/dcd format. 
Example:
```
python automatization_MDpocket -t [path_to_trajectory_filename] -p [path_to_pdb_file]
```

#### OPTIONAL ARGUMENTS

##### -i, --isovalue:
Isovalue selected for the extraction of pocket coordinates. Higher isovalues will have pockets that are more maintained during the trajectory. Otherwise, lower isovalues will select transient pockets. Default: isovalue = 3. Using an isovalue of 3 we are able to obtain transient pockets.

##### -s, --pdb_step:
The trajectory is reduced in order to decrease the computational cost of the algorithm. In order to do so it is only used 1/pdb-step frames. So, the trajectory is reduced by 1/pdb-step. Default: pdb_step = 5. Using a pdb_step = 5, we are able to reduce considerably the computation time without affecting the results obtained by MDpocket.

##### -c, --cpu:
Number of CPUs used when computing the descriptors of each pocket. Using more CPUs will reduce the computational time of the algorithm. Default = it will use all the CPUs of your computer. If you don't want to use all of them you can specify the number of cpus as 3 using:
```console
python automatization_MDpocket -d [path_to_directory] -c 3
```

### WORKFLOW

* Matching pdb files with trajectories: First of all as mentioned before the algorithm matches each trajectory with its corresponding model file (pdb).

* Creating MDpocket input: MDpocket uses a list of pdb files as input. For this reason each frame of the trajectory is converted into a pdb file. We will create a folder named snapshots_input_NameTrajectoryFile where all pdb files will be stored
  with the name: md + number of the snapshot in the simulation + .pdb. For a matter of computation time it will just be created a 1 pdb file every 5 frames but it can be changed with the variable pdb_step.
  The pdb files created will only contain the protein. Finally, it will be created a file (snapshot_list_file_TrajectoryFileName) where the path to every pdb file is found. This file will be used as input for MDpocket.
  
* Running MDpocket: MDpocket is run over the simulation and its output is stored in the folder mdpocket. To see more information about the output of MDpocket you can go [here](http://fpocket.sourceforge.net/manual_fpocket2.pdf) (page 21-24).

* Extracting pocket coordinates: from the output file of MDpocket mdpout_freq.dx which has a grid that contains a measure of frequency of how many times the pocket was open during a MD trajectory. The frequency of the pocket is expressed as the number of Voroni Vertices in a 8 Armstrong cube around each grid point per snapshot, so the more a cavity is conserved the highers is the value. In order to detect transient pockets it is extracted the grid points with an isovalue equal or higher than 3. As a result, we obtain all the coordinates detected as a pockets with an isovalue equal or higher than 3. These coordinates can be visualized with the following file: Extract_ISO_output_3_TrajectoryFileName.
  
* Clustering coordinates into pockets: Using density-based spatial clustering of applications with noise (DBSCAN) algorithm the coordinates are clustered into different pockets. The idea behind DBSCAN is that a point belongs to a cluster if it is close to many points from that cluster. In addition, it is albe to detect outliers (coordinates that are far from the others) allowing not to take into consideration coordinate points that are outliers. The set of coordinate points that are creating the pocket are stored in a pdb file. The different pdb files with the coordinates of each pocket are stored in DBSCANclustering_coordinates_pdb_TrajecoryFileName.
* 2nd run of MDpocket: The pdb files obtained in the previous step will be used to run for a second time MDpocket in order to obtain the descriptors of each pocket (Volume, hydrophobicity score...) that will be stored in the descriptorPocekts_TrajectoryFileName. In addition, it will be create a graphic representing the Volume of the pocket along the trajectory and they will be stored in Volume_graphic_TrajectoryFileName.

## OUTPUT

* DBSCANclustering_coordinates_pdb_TrajectoryName. This directory contains pdb files of the different pockets found in the trajectory. If you want to load all the pockets on VMD just open the TKconsole and type this command on the DBSCANclustering_coordinates directory:
```
set pdblist [glob *pdb]
foreach pdb $pdblist {
mol new $pdb
}
```

If you also open the trajectory, this will allow you to see all the pockets found and see which pocket is the one that you are more interested in.


* snapshots_TrajectoryName: in this folder will be stored the pdb files of each frame of the trajectory as it is the required input format for MDpocket.
* mdpocket: Here are stored the files created by the first run of MDpocket.  
    * mdpout_dens.dx: A grid is superimposed to all alpha spheres of all snapshots and the number of alpha spheres around each grid point is counted. In order to see more clearly it is recommended to change the representation (using VMD going to Graphics -> Representations). Here you can change the isovalue (higher more conserved pockets). Here an image of how to change representations to see pockets on MD trajectory usng VMD:
    ![VMD conf](https://user-images.githubusercontent.com/57498211/161777199-65b3a424-e69f-4e31-b19c-2afd30fa8803.png) 
    * mdpout_freq.dx: This file is very similar to the previous grid file. Here the grid contains a measure of frequency of how many times the pocket was open during a MD trajectory averaged by the number of snapshots. Thus, this gives a range of possible iso-values between 0-1. 
    * mdpout_dens_iso_8.pdb: This file contains all grid points having 3 or more Voronoi Vertices in the 8A3 volume around the grid  point for each snapshot. 
    * mdpout_freq_iso_0_5.pdb: This is similar to the previous pdb file, just being produced on the frequency grid with a cutoff of 0.5.
    




* MDpocket_Voroni_Vertices_TrajectoryName: In this directory you can find the trajectory file of the pockets during the trajectory. Using these files you will be able to see the movement of the pocket. Be careful, VMD does not read this file, as from one snapshot to the other a different number and type of Voronoi vertices can be part of the model.
* MDpocket_atoms_TrajectoryName: Here, you can find pdb files containing all receptor atoms defining the binding pocket in each frame analyzed. 
* descriptorPockets_TrajectoryName: Last but not least, here you will find the descriptors of each pocket (Volume, mean local hydrophobic density, mean alpha sphere solvent accessibility...). 
This output file can be easily analyzed using R, gnuplot or other suitable software. An example R output for the pocket volume would be like:
```
r=read.table("pocket_num_0_descriptors.txt",h=T)
ylim=c(400,1200)
plot(r[,"pock_volume"],ty='l',ylim=ylim,main="",xlab="",ylab="")
par(new=T)
plot(smooth.spline(r[,"pock_volume"],df=40),col="red",lwd=3,ylim=ylim,ty
="l",xlab="snapshot",ylab="volume")
```
With this, you should be able to see if the volume of the pocket increases or decreases during the trajectory.

* Volume_graphic_TrajectoryName: Here there is a graphical representation of the volume of each pocket along the trajectory in order to easily see the opening and clousure of the pocket if it is required.
* Final_OUTPUT_TrajectoryName: Finally a txt file is created with the most important features of each pocket in order to easily detect the most interesting pockets detected over the trajectory. In this file we can see the volume, average volume,  atoms, the mean local hydrophobic density, average Polarity score and the average hydrophobicity score which are the descriptors used for MDpocket to calculate the druggability score of each pocket.

## ADDITIONAL INFORMATION

* DRUGGABILITY SCORE:
To obtain a druggability score in MDpocket, you need to run MDpocket using the -S flag. However, using the -S flag, less pockets are detected and consequently we are losing information about transient pockets. For this reason the -S flag is not used in this algorithm and consequently we are not geting the druggabilit score. 
However, the formula used to obtain the druggability score is the following one (obtained from [here](https://www.researchgate.net/publication/45504065_Understanding_and_Predicting_Druggability_A_High-Throughput_Method_for_Detection_of_Drug_Binding_Sites): 


![Druggability score formula](https://user-images.githubusercontent.com/57498211/161771891-dc589b56-2b9e-4658-9df1-343995802bd7.png)
![Descriptors_druggability_score](https://user-images.githubusercontent.com/57498211/161772411-d3b6b0fa-2026-4f48-ac2f-b621f4274825.png)

However this formula can't be applied to the pocket descriptors as the characterization step is delimiting the pocket by a user input, calculating a druggabiity score isn't relevant as the score hasn't been trained on such input, it is thus deactivated.
However, an alternative option given by the creators of MDpocket is to use the mean local hydrophobic density as an aproximation of the druggability score as it is the descriptor that contributes more to the druggability score. The mean local hydrophobic density of the binding site is considered as the most predictive descriptor as it combines the size and spatial distrubution of hydrophobic agglomerations into a single number So, the mean local hydrophobic density is a predictive descriptor of interest in order to know the druggability of the pocket (the higher the more druggable).


