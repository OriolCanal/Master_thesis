#from datetime import datetime
import os, os.path
import re
import collections
import math
import fnmatch
from subprocess import run, PIPE, Popen
from scipy.interpolate import make_interp_spline, BSpline
import pathlib
import json
import subprocess
import numpy as np
import matplotlib.pyplot as plt
#%matplotlib inline
from mpl_toolkits.mplot3d import Axes3D
from sklearn.cluster import DBSCAN
from sklearn import preprocessing
import sys
import multiprocessing
from functools import partial
import argparse

# Bokeh libraries
from bokeh.io import output_file, output_notebook
from bokeh.plotting import figure, show, save
from bokeh.models import ColumnDataSource
from bokeh.layouts import row, column, gridplot
from bokeh.models.widgets import Tabs, Panel

#begin_time = datetime.now()
parser = argparse.ArgumentParser(description = "Program to autmatically detect pockets and its decriptors over MD simulations using MDpocket.")
group = parser.add_mutually_exclusive_group(required=True)
group.add_argument ("-d", "--input_directory", action = "store" ,help="directory where the pdb files and trajectories are stored")
group.add_argument ("-t", "--trajectory_file", action = "store", help ="trajectory filename (xtc or dcd formats)" )
parser.add_argument ("-p", "--pdb_file", action = "store", required = False, help = "model pdb file of the protein") 
parser.add_argument ("-i", "--isovalue", action = "store", required = False, default = 2.5, type = float, help = "MDpocket isovalue to extract pockets (higher isovalue will consider more conserved pockets). Recommended isovalue between 2 (transient) and 5 (more conserved).")
parser.add_argument ("-s", "--pdb_step", action = "store", required = False, default = "5", type = str , help = "Reduction of the molecular dynamics frames to decrease the computation time" )
parser.add_argument ("-c", "--cpu", action = "store", required = False, default = 0, type = int, help = "Number of cpus to use in the parallization step.")
parser.add_argument ("-a", "--structural_alignment", action = "store", required = False, choices = ['yes', 'no'], default = 'no', help ="if the trajectory is not previously aligned, include -a yes")
args = parser.parse_args()



input_directory = args.input_directory
print (input_directory)
isovalue = args.isovalue
pdb_step = args.pdb_step
cpu = args.cpu
alignment = args.structural_alignment
input_trajectory_file= args.trajectory_file
input_topology_filename = args.pdb_file

def matching_pdb_traj(directory):

  """This function matches de model file (pdb) with its trajectory in a directory."""
  directory = str(directory)
  trajectories_analysed = []
  #Find the pdb files in the  directory
  pdb_files = [f for f in os.listdir(directory) if f.endswith(".pdb")]
  #Find trajectory files ("xtc" or "dcd") in the directory
  trajectory_files = [f for f in os.listdir(directory) if f.endswith(".xtc") or f.endswith(".dcd")]
  trajectories_folder_analysed =  [f for f in os.listdir(directory) if f.endswith("_mdpocket")]
  for trajectory_folder in trajectories_folder_analysed:
    trajectory = trajectory_folder.replace("_mdpocket","")
    trajectories_analysed.append(trajectory)
    
  for trajectory_analysed in trajectories_analysed:
    trajectory_files.remove(trajectory_analysed)
  #print(pdb_files)
  #dictionary where pdb will be matched with the trajectory
  id_dic = {}
  #print (trajectory_files)
  for f in pdb_files:
    if os.path.isfile(f):
      pdbid_regex=re.compile(r"dyn_[0-9]+")
      pdb_id = pdbid_regex.search(f)
      #print (pdb_id[0])
      pdb_num_regex = re.compile(r"[0-9]+")
      pdb_num = pdb_num_regex.search(pdb_id[0])
      #print(pdb_num[0])
      if f not in id_dic:
          id_dic[f] = list()
    for g in trajectory_files:
      if os.path.isfile(g):
        #print( "trajectory file:", g)
        trajid_regex = re.compile(r"trj_[0-9]+")
        traj_id = trajid_regex.search(g)
        traj_num_regex = re.compile (r"[0-9]+")
        traj_num = traj_num_regex.search(traj_id[0])
        #print (traj_num[0])

      if pdb_num[0] == traj_num[0]:
          id_dic[f].append(g)
      
  #print(id_dic)
  return(id_dic)



#Perform the pockets analysis
def extract_ISO (
  inputfile : str,
  pathOutput : str,
  iso_value : float):
  

  f=open(inputfile,"r")


#get the axis that shows the most variation during the trajectory, this will be the leading axis

#read the header - here is an example
  header=""
  tmp=f.readline()
  while tmp[0]!="o" : 
    header= header + tmp
    tmp=f.readline()
#print header

#read the grid size
  r=re.compile('\w+')
  gsize=r.findall(tmp)
  gsize=[int(gsize[-3]),int(gsize[-2]),int(gsize[-1])]
#print gsize

#read the origin of the system
  line=f.readline().split()
  origin=[float(line[-3]),float(line[-2]),float(line[-1])]
#print origin

#read grid space
  line=f.readline().split()
  deltax=[float(line[-3]),float(line[-2]),float(line[-1])]
  line=f.readline().split()
  deltay=[float(line[-3]),float(line[-2]),float(line[-1])]
  line=f.readline().split()
  deltaz=[float(line[-3]),float(line[-2]),float(line[-1])]


#pay attention here, this assumes always orthogonal normalized space, but normally it should be ok
  delta=np.array([deltax[0],deltay[1],deltaz[2]])

#read the number of data
  f.readline()
  r=re.compile('\d+')
  n_entries=int(r.findall(f.readline())[2])

  if(n_entries!=gsize[0]*gsize[1]*gsize[2]) : sys.exit("Error reading the file. The number of expected data points does not correspond to the number of labeled data points in the header.")

#create a 3D numpy array filled up with 0


#initiate xyz counter for reading the grid data
  z=0
  y=0
  x=0

  print("Reading the grid. Depending on the number of data points you have this might take a while....")
  path=open(pathOutput,"w")

  counter=1
  for count in range(n_entries//3) :
    c=f.readline().split()
    if(len(c)!=3) : 
      print("error reading grid data")
      sys.exit("exiting the program")
    for i in range(3):
      if (iso_value<0 and float(c[i]) < iso_value) or (iso_value > 0 and float(c[i]) > iso_value) :
        path.write('ATOM  %5d  C   PTH     1    %8.3f%8.3f%8.3f%6.2f%6.2f\n'%(counter,origin[0]+float(x)*delta[0],origin[1]+float(y)*delta[1],origin[2]+float(z)*delta[2],0.0,0.0))
        counter+=1
      z+=1
      if z >= gsize[2]:
        z=0
        y+=1
        if y >=gsize[1]:
          y=0
          x+=1


  path.close()
  f.close()

  print("finished writing %s"%(pathOutput))


######### END OF EXTRACT_ISO FUNCTION ##############

  




##############CREATE MD POCKET INPUT FUNCTION #############################

def creteMDpocket_input(
  input_trajectory_filename : str,
  input_topology_filename : str,
  pdb_step : str,
  outputFile: str):
  """ It creates a pdb file for every snapshot of the trajectory given as input.
  We need to give to the function the topology file and the a topology file and 
  the step between frames of the new reduced trajectory.
  We will create a folder named snapshots_input_trajectory_filename where all pdb files will be stored
  with the name: md + number of the snapshot in the simulation + .pdb 
  The pdb files created will only contain the protein"""


  if alignment == 'yes':
    # name of the trajectory file aligned
    #### UNCOMMENT THE 9 LINES BELOW IF STRUCTURAL ALIGNMENT OF THE TRAJECTORY IS NEEDED. 
    input_trajectory_aligned = "aligned_" + input_trajectory_filename 
  
  
    print("PERFORMING THE STRUCTURAL ALIGNMENT OF THE TRAJECTORY")
    #first let's do a structural alignment. This step is necessary (reasons indicated on the MDpocket paper.) 
    p = subprocess.Popen(['gmx', 'trjconv',
                          '-s', input_topology_filename, '-f', input_trajectory_filename,
                          '-o', input_trajectory_aligned, '-fit', 'rot+trans'],
                          stdin=subprocess.PIPE)
    p.communicate(b'3\n0\n') # alignment on c-alpha, output everything
    p.wait()
    input_trajectory_filename = input_trajectory_aligned



  #Creating a folder where all pdb snapshots will be stored
  snapshots_folder = "snapshots_" + input_trajectory_filename
  if not os.path.exists(snapshots_folder):
    logs = run([
          "mkdir",
          snapshots_folder,
          ], stdout=PIPE).stdout.decode()

    
    
  
    #Obtain pdbs from the trajectory using gmx trjconv: -s (topology), -f (trajectory), -o (output pdb name) -nzero -skit(to indicate pdb-steps over the simulation)

  snapshots_output = snapshots_folder + '/md.pdb'

  
  print("Creating pdb snapshots from trajectory file, they will be stored in:" + snapshots_folder)
  

  #if the trajectory is a dcd file it will be converted to a xtc file in order to perform the analysis.
  if input_trajectory_filename.endswith("dcd"):
    input_trajectory_xtc = input_trajectory_filename.replace("dcd","xtc")

    p = subprocess.Popen(['mdconvert', input_trajectory_filename, "-o", input_trajectory_xtc],
                       stdin=subprocess.PIPE)
    p.wait()
    input_trajectory_filename = input_trajectory_xtc


#command to obtain pdb files from the trajectory. The echo 1 is to indicate that we only want to store the protein in the pdb files.
    # gmx trjconv -s " + input_topology_filename + " -f " + input_trajectory_aligned + " -sep -nzero 1 -skip " + pdb_step + " -o " + snapshots_output
  p = subprocess.Popen(['gmx', 'trjconv',
                        '-s', input_topology_filename, '-f', input_trajectory_filename,
                        '-sep', '-nzero', '1', '-skip', pdb_step, '-o', snapshots_output],
                       stdin=subprocess.PIPE)
  p.communicate(b'Protein\n') # We only obtain pdb snapshots of the protein
  p.wait()

  




  #just a helper function to do some stupid concatenation
  def getFname(s):
    return (snapshots_folder+os.sep+s)

  #get the fpocket output folders
  snapshots=fnmatch.filter(os.listdir(snapshots_folder),"*.pdb") 

  snapshots=[os.path.abspath(getFname(sn)) for sn in snapshots] 

  #set a pattern to find digits
  RE_DIGIT = re.compile(r'(\d+)')   
  #create on the fly function (lambda) to return numbers in filename strings 
  ALPHANUM_KEY = lambda s: [int(g) if g.isdigit() else g for g in RE_DIGIT.split(s)]
        
  snapshots.sort(key=ALPHANUM_KEY)      #sort by these numbers in filenames

  fout=open(outputFile,"w")
  [fout.write(sn+"\n") for sn in snapshots]
  fout.close()

  number_snapshots = len([name for name in os.listdir('./' + snapshots_folder) if os.path.isfile(name)])
  
  




####################CREATE MDPOCKET INPUT FUNCTION UNTIL HERE#########################


###############FUNCTION TO RUN MDPOCKET TO GET THE CHARACTERISTICS OF THE POCKETS (MDPOCKET 2ND ROUND)#############
#defining the variable pockets trajectory to be used in the function


def run_mdpocket_2nd(pockets_trajectory, pocketname):
  global new_pdb_filename
  global outputFile
  new_pdb_filename =  str(pocketname) + '_coordinates_DBSCAN.pdb'
  outputFile= "snapshot_list_file_"  +pockets_trajectory + ".txt"
  #Finally we run MDpocket for a second time using each pocketID.pdb created to obtain the descriptors of each pocket.
  cmd_pocket_descriptors = "mdpocket -L " + outputFile + " -o " +"pocket_num_" + str(pocketname) + " -f " + new_pdb_filename + " -v 10000"
  print("ANALYSING POCKET " + str(pocketname))
  os.system(cmd_pocket_descriptors)




#######################APPLYING DBSCAN FUNCTION ##########################



def applying_DBSCAN(
  grid_filename: str,
  output_analysis_filename : str,
  input_trajectory_filename : str,
  pdb_step: int):

  """ It takes as input the pdb of the coordinates that have an isovalue over the cutoff and 
  applies DBSCAN to it in order to cluster the coordinates into different pockets."""
  

  output_analysis = []
  coordinates_index= 0


  # Lista -> will be stored the x,y,z coordinates of each point
  lista = []

  # Atom_num -> will store the atom number of each coordinate. Uri: I think we don't need that
  Atom_num = []
  

  #reads the coordinates file and saves each coordinate as a float
  with open(grid_filename) as f:
    lines = f.readlines()

  for line in lines:
    x_coordinates=float(line[31:38].strip())
    y_coordinates=float(line[39:46].strip())
    z_coordinates=float(line[47:54].strip())
    ATOM_number = line[7:11].strip()


    Atom_num.append(ATOM_number.strip())       


    coordinates_index = coordinates_index +1


    lista.append([x_coordinates, y_coordinates, z_coordinates])

  # the coordinates will be saves as a numpy array in data variable as is the input needed to apply DBSCAN. 
  data = np.array(lista)
  #print (data)

  #plot the coordinates found as pocket in a 3D plot.(commented, 
  #if necessary it can be produced as output)

  #fig = plt.figure()
  #ax = Axes3D(fig)
  #ax.scatter(data[:,0], data[:,1], data[:,2], s=300)
  #ax.view_init(azim=200)
  #plt.show()


  #apply DBSCAN to the data. Here the important parameters that can be changed are the eps and the min samples
  model = DBSCAN(eps=1.5, min_samples=8)
  model.fit_predict(data)
  pred = model.fit_predict(data)


  #plotting a 3D plot where each coordinated found as pocket is coloured depending on the cluster they belong.
  #fig = plt.figure()
  #ax = Axes3D(fig)
  #ax.scatter(data[:,0], data[:,1], data[:,2], c=model.labels_, s=300)
  #ax.view_init(azim=200)
  #plt.show()
  


  i=0

  #Applying DBSCAN we get each coordinate grouped in model.labels that we will use to know which coordinates are clustering to one pocket 
  pocketID = model.labels_



  #Obtaining the number of pockets found using DBSCAN
  num_pocketID = max(pocketID)


  if num_pocketID == -1 or num_pocketID == 0:
    print('WARNING: No pockets were found')
  else:
    print("The total number of pockets found is: " + str(num_pocketID))


  #In coordinate_pocketID we will store the coordinates and the pocketID for each coordinate
  coordinate_pocketID = []
  for coordinate in lista:
  
    coordinate_pocketID.append((coordinate, pocketID[i]))
    i = i+1
  #print (coordinate_pocketID)

  #We loop onto the different pocketsIDs and we take all the coordinates that have the same 
  #pocketID and it is created a pdb file from the coordinates that are clustered in the same pocketID
  # It doesn't consider the outliers because the outliers have the pocketname -1
  pocketname = 0
  for pocketname in range(num_pocketID+1):
    new_pdb_filename =  str(pocketname) + '_coordinates_DBSCAN.pdb'
    new_pdb_lines= []
    pocket_coordinates=[]
    lines_count=0
    for coordinates, pocketIDs  in coordinate_pocketID:
      if pocketIDs == pocketname:
        pocket_cooridnates=pocket_coordinates.append(coordinates)
    for x,y,z in pocket_coordinates:
      #We save the coordnates in this format to crete the PDB file with a correct format
      x_cord = str(x).rjust(8, ' ')
      y_coord = str(y).rjust(8,' ')
      z_coord = str(z).rjust(8,' ')
      lines_count += 1
      atom_num = str(lines_count).rjust(6,' ')
      line = "ATOM "+ atom_num +"  C   PTH     1    "+ x_cord + y_coord + z_coord + "  0.00  0.00\n"  
      new_pdb_lines.append(line)


    # The pdb file created will be named pocketID.pdb (e.g. 1.pdb will have the coordinates clustered in the pocke number 1)
    with open(new_pdb_filename,'w') as file:
      for line in new_pdb_lines:
        file.write(line)
    outputFile= "snapshot_list_file_" + input_trajectory_filename + ".txt"
  

  #paralelization of the run_mdpocket_2nd function to make the slowest step of the algorithm faster. 
  iterable = range(0,int(num_pocketID)+1)
  if cpu == 0:
    pool_obj = multiprocessing.Pool(6)
  else:
    pool_obj = multiprocessing.Pool(cpu)

  # to use pool.map() with more than one parameter is used functools.partial
  func = partial(run_mdpocket_2nd, input_trajectory_filename)
  pool_obj.map(func, iterable)
  pool_obj.close()
  pool_obj.join()


    #Mine data from the mdpocket descriptors output file
  
  #Variables where pocket descriptors will be stored to be able to calculate the druggability score.
  average_polarity_score_list = []
  average_hydrophobicityDensity_list = []
  average_hydrophobicityScore_list = []
  drugscore_list = []

  for pocketname in range(num_pocketID+1):

    descriptors_data = {}
    mdpocket_descriptors_filename = "pocket_num_" + str(pocketname) + "_descriptors.txt"
    with open(mdpocket_descriptors_filename,'r') as file:
      # The '[:-1]' is to remove the break line at the end of each line
      entries = re.split("[ ]+", next(file)[:-1])
      for entry in entries:
        descriptors_data[entry] = []
      for line in file:
        line_data = re.split("[ ]+", line[:-1])
        for i, value in enumerate(line_data):
          descriptors_data[entries[i]].append(value)

    # Mine the atoms implicated in each pocket
    # In this file atoms are listed for each frame, but they are always the same
    # For this reason, we mine only atoms in the first frame
    atoms = []
    atoms_filename = "pocket_num_" + str(pocketname) + "_mdpocket_atoms.pdb"
    with open(atoms_filename,'r') as file:
      for line in file:
        line_data = re.split("[ ]+", line[:-1])
        if line_data[0] == 'MODEL':
          continue
        if line_data[0] == 'ATOM':
          atoms.append(int(line_data[1]))
        if line_data[0] == 'ENDMDL':
          break

    # Format the mined data and append it to the output data
    # NEVER FORGET: The 'descriptors_data' object contains a lot of data, not only pocket volumes
    # (e.g. drugability score)
    
      #we define the snapshos folder in this function again because it wasn't defined in this function before
    snapshots_folder = "snapshots_" + input_trajectory_filename
    #counting the number of files are created to know the number of frames we analysed:
    count = 0
    for path in pathlib.Path("./"+snapshots_folder).iterdir():
      if path.is_file():
        count += 1

   #PLOTTING y = VOLUME OF EACH POCKET VS FRAME
    frames = list(range(0,count))
    real_frames = [element * int(pdb_step) for element in frames]
    volume = list(map(float, descriptors_data['pock_volume']))

    # The pocket will be considered as transient if in 10 consecutive frames the volume of the pocket = 0.
    transient = False
    for i,j,q,k,l,m,n,o,p,q in zip(volume, volume[1:],volume[2:],volume[3:],volume[4:],volume[5:],volume[6:],volume[7:],volume[8:],volume[9:]):

      #checking if 10 consecutive volumes = 0
      if (i + j + q + k + l+m+n+o+p+q) == 0:
        transient = True
        break #If 10 consecutive 0 are found, we don't need to continue the loop, the pocket will be transient.

    #Average volume of the pocket over the frames to add in the plot if it is needed.
    average_volume = round(sum(volume)/len(volume),2)

    # Determine where the visualization will be rendered
    output_file('volume_pocket_'+str(pocketname)+'.html')# Render to static HTML, or
    #output_notebook() # Render inline in a Jupyter Notebook (uncomment it if necessary)

    #Set up the figures
    fig = figure(title= 'Pocket ' + str(pocketname) + ' volume',
    x_axis_label = 'Frame', y_axis_label = 'Volume')
    fig.line(x=real_frames, y= volume,
      color = 'red', line_width=1,
      legend = 'Pocket ' +str(pocketname) + ' avg volume = ' + str(average_volume))

    #figure legend object on top left of the graph:
    fig.legend.location = 'top_left'

    save(fig)
    # Saving the variables that we need to calculate the drugability score
    # The coefficients will be stored as they have to be normalized after

    hydrophobicity_density = list(map(float, descriptors_data['mean_loc_hyd_dens']))
    hydrophobicity_density = [i for i in hydrophobicity_density if i != 0]
    if len(hydrophobicity_density) != 0:
      average_hydrophobicityDensity = round(sum(hydrophobicity_density)/len(hydrophobicity_density),2)
    else:
      average_hydrophobicityDensity = 0
    #np_average_hydrophobicityDensity = np.array(average_hydrophobicityDensity_list).reshape(-1,1)


    hydrophobicity_score = list(map(float, descriptors_data['hydrophobicity_score']))
    hydrophobicity_score = [i for i in hydrophobicity_score if i != 0]
    if len(hydrophobicity_score) != 0:
      average_hydrophobicityScore = round(sum(hydrophobicity_score)/len(hydrophobicity_score),2)
    else:
      average_hydrophobicityScore = 0


    polarity_score = list(map(float, descriptors_data['polarity_score'])) 
    polarity_score = [i for i in polarity_score if i != 0]
    if len(polarity_score) != 0:
      average_polarity_score = round(sum(polarity_score)/len(polarity_score),2)
    #np_average_polarity_score = np.array(average_polarity_score_list).reshape(-1,1)
    else:
      average_polarity_score = 0

    output = {
      'name': pocketname,
      'volumes': list(map(float, descriptors_data['pock_volume'])),
      'average_volume': average_volume,
      'atoms': atoms,
      'average_Polarity_score' : average_polarity_score,
      'average_hydrophobicity_score' : average_hydrophobicityScore,
      'average_hydrophobicity_Density': average_hydrophobicityDensity,
      'transient' : transient,
    }


    #Set the dict where all ouptut data will be stored and append the output previously defined.
    output_analysis.append(output)



    # Export the analysis in json format
    with open(output_analysis_filename, 'w') as file:
        json.dump({ 'data': output_analysis}, file)

  #not useful but I let it here in case is necessary in future
  #PREDICTING DRUGGABILITY SCORE
    #The druggability score will be predicted after using the scoring function of mdpocket that can be observed in the following paper:
    # https://www.researchgate.net/publication/45504065_Understanding_and_Predicting_Druggability_A_High-Throughput_Method_for_Detection_of_Drug_Binding_Sites
    # The formulas can be found in page 9 of the article. The coefficient ss found in table 3 are the beta coefficients used in the formula.
  
  #scaler = preprocessing.MinMaxScaler()
  #normalized_hydrophobicityDensity = scaler.fit_transform(np_average_hydrophobicityDensity)
  #normalized_polarity_score = scaler.fit_transform(np_average_polarity_score)
  #print("normalized polarity score = ", normalized_polarity_score, "normalized hydrophobicity density = ", normalized_hydrophobicityDensity,
  #  "average hydrophobicity score = ", average_hydrophobicityScore_list)
    
  #iterating over the number of pockets:
  #for x in range(num_pocketID):
  #  print ("x = ", x)
  #  f1d1 = (math.e**(-(-5.141)+6.579*float(normalized_hydrophobicityDensity[x])))/(1+math.e**(-(-5.141)+6.579*float(normalized_hydrophobicityDensity[x])))
  #  f2d2 = (math.e**(-(-2.669)+0.056*float(average_hydrophobicityScore_list[x])))/(1+math.e**(-(-2.669)+0.056*float(average_hydrophobicityScore_list[x])))
  #  f3d3 = (math.e**(-(-2.445)+2.762*float(normalized_polarity_score[x])))/(1+math.e**(-(-2.445)+2.762*float(normalized_polarity_score[x])))
  #  print("f1d1", f1d1)
  #  print("f2d2", f2d2)
  #  print("f3d3", f3d3)
    
    #z = -6.238 + 4.592 * f1d1 + 5.717 * f2d2 + 3.985 * f3d3

  #  print("z:", z)

  #  drugscore = (math.e**(-(float(z)))/(1+math.e**(-(float(z)))))
  #  print("drugscore:", drugscore)
  #  drugscore_list.append(drugscore)
  #  print ("drugscore list", drugscore_list)

    #with open(output_analysis_filename, 'a') as file:
      #json.dump({ 'druggability_score_pocket_'+ str(x): drugscore}, file)
#################################################################################



##########################MAIN FUNCTION TO RUN MDPOCKET##################################

def pockets (
  input_topology_filename : str,
  input_trajectory_filename: str,
  #output_analysis_filename : str,  ###remove the # when we go for the output step
  isovalue: float,
  pdb_step: str):

  #Defining the file where the list of pdb files will be stored
  outputFile = "snapshot_list_file_" + input_trajectory_filename + ".txt"

  print (outputFile)
  creteMDpocket_input(
  input_trajectory_filename,
  input_topology_filename,
  pdb_step,
  outputFile)



  

  #Create a new folder to store all output files so they do not overcrowd the main directory

  mdpocket_folder = 'mdpocket'
  if not os.path.exists(mdpocket_folder):
    logs = run([
      "mkdir",
      mdpocket_folder,
      ], stdout=PIPE).stdout.decode()

  print ("mdpocket_folder: " , mdpocket_folder)
  #Run MDpocket to find new pockets over the simulation

  mdpocket_output = mdpocket_folder + '/mdpout'

  print ("mdpocket_output: " , mdpocket_output)

  #Set the filename of the mdpocket output we are interested in

  grid_filename = mdpocket_output + '_dens.dx'

  #Skip the last step if  the output file already exists and is not empty
  if not os.path.exists(grid_filename) or os.path.getsize(grid_filename) == 0:
    mdpocket_run = "mdpocket -L " + outputFile + " -o " + mdpocket_output
    os.system(mdpocket_run)

    #-----------------------------------------------
#We use extract_ISO function to extract the coordinates that have isovalue over the cutoff isovalue selected
  extract_ISO (
  grid_filename,
  "Extract_ISO_output_" + str(isovalue) + "_" + input_trajectory_filename ,
  isovalue)


  applying_DBSCAN(
  "Extract_ISO_output_" + str(isovalue) + "_" + input_trajectory_filename,
  "Final_OUTPUT_" + input_trajectory_filename,
  input_trajectory_filename,
  pdb_step)


  #Creting folders to store output:

  atom_folder = 'MDpocket_atoms_' + input_trajectory_filename 
  if not os.path.exists(atom_folder):
    logs = run([
      "mkdir",
      atom_folder,
      ], stdout=PIPE).stdout.decode()

      # In the atom_folder will be stored a pdb file containing all receptor atoms surrouding the
      # selected pocket region. It is a NMR like file containing each snapshot as seperated model.
      # Can be viewed using VMD and PyMOL.
    cmd_atoms = "mv pocket_num_*_mdpocket_atoms.pdb ./" + atom_folder
    os.system(cmd_atoms)


  #In the pocket_descriptor_folder will be stored a txt file containing the fpocket descrptors
  #of the selected pocktes region for each snapshot. To analyse using R.
  pocket_descriptor_folder = 'descriptorPockets_' + input_trajectory_filename 
  if not os.path.exists(pocket_descriptor_folder):
    logs = run([
      "mkdir",
      pocket_descriptor_folder,
      ], stdout=PIPE).stdout.decode()

    cmd_descriptors = "mv *descriptors.txt ./" + pocket_descriptor_folder
    os.system(cmd_descriptors)



  #In the Voroni_vertices_folder will be stored te pdb files contaning all Vooni vertices within
  #the selected pocket region for all snapshots. This file s a NMR like file, containing each
  #snapshots as separated model. Can be viewed using PyMOL.
  pocket_Voroni_vertices_folder = 'MDpocket_Voroni_Vertices_' + input_trajectory_filename 
  if not os.path.exists(pocket_Voroni_vertices_folder):
    logs = run([
      "mkdir",
      pocket_Voroni_vertices_folder,
      ], stdout=PIPE).stdout.decode()

    cmd_Voroni_Vertices = "mv pocket_num_*_mdpocket.pdb ./" + pocket_Voroni_vertices_folder
    os.system(cmd_Voroni_Vertices)



  #In the cluster_pockets_DBSCAN_folder will be stored the pdb files with the coordinates of each pocket that are obtained as output 
  #after applying DBSCAN_clustering. 
  cluster_pockets_DBSCAN_folder = 'DBSCANclustering_coordinates_pdb_' + input_trajectory_filename 
  if not os.path.exists(cluster_pockets_DBSCAN_folder):
    logs = run([
      "mkdir",
      cluster_pockets_DBSCAN_folder,
      ], stdout=PIPE).stdout.decode()

    cmd_pockets_DBSCAN = "mv *_coordinates_DBSCAN.pdb ./" + cluster_pockets_DBSCAN_folder
    os.system(cmd_pockets_DBSCAN)

  #In the volume_graphics folder it will be stored the graphic of the volumes of the pockets over the frames
  volume_graphic_folder = 'Volume_graphic_' + input_trajectory_filename
  if not os.path.exists(volume_graphic_folder):
    logs = run([
      "mkdir",
      volume_graphic_folder,
      ], stdout=PIPE).stdout.decode()

    cmd_graphic_folder = "mv volume_pocket_* ./" + volume_graphic_folder
    os.system(cmd_graphic_folder)


  #Creating a directory where all the folders and files created until now will be stored.
  general_directory =  input_trajectory_filename + '_mdpocket'
  if not os.path.exists(general_directory):
    logs = run([
      "mkdir",
      general_directory,
      ], stdout=PIPE).stdout.decode()


  #MOVING ALL THE FILES AND FOLDERS GENERATED DURING THE MDPOCKET RUN INTO THE GENERAL DIRECTORY
  cmd_general_directory = "mv *" + input_trajectory_filename + "* ./" + general_directory
  os.system(cmd_general_directory)
  cmd_mv_mdpocket_General = "mv mdpocket ./" + general_directory
  os.system(cmd_mv_mdpocket_General)
  cmd_traj_to_initial_place = "mv ./" + general_directory + "/" + input_trajectory_filename + " ."
  os.system(cmd_traj_to_initial_place)
  #print(datetime.now() - begin_time)


def get_headlines_gridfile (grid_filename = str):
  """It returns the header lines of the grid file. In any case it is needed the origin and dimensions of the system."""
  grid_values_pattern = "^([.0-9]+) ([.0-9]+) ([.0-9]+) $"
  # First of all, get all header lines from the original grid file
  # We need them to write grid files further
  with open(grid_filename,'r') as file:
    header_lines = []
    header_pattern = "^[a-z]"
    for line in file:
      if re.match(header_pattern, line):
        header_lines.append(line)
      elif re.match(grid_values_pattern, line):
        break
  return(header_lines)



if __name__ == "__main__":
  if input_directory != None:
    dictionary = matching_pdb_traj(input_directory)
    for model, trajectories in dictionary.items():
      for trajectory in trajectories:
        pockets(
          str(model),
          str(trajectory),
          isovalue,
          pdb_step)

  else:
    pockets(
          str(input_topology_filename),
          str(input_trajectory_file),
          isovalue,
          pdb_step)





























