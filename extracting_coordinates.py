with open("extractIso_output_isovalue_5.pdb ") as f:
    lines = f.readlines()
charged_res_coord = [] # store x,y,z of extracted charged resiudes 
for line in lines:
    if line.startswith('ATOM'):     
    atom_coord.append(line)

  for i in range(len(atom_coord)):
   for item in charged_res:
     if item in atom_coord[i]:
        charged_res_coord.append(atom_coord[i].split()[1:9])

print (charged_res_coord)