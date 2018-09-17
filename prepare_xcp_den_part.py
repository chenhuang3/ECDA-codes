#!/usr/bin/python
import os 
import numpy as np

os.system("rm -rf subsys?")
os.system("rm -fr subsys??")
os.system("rm -fr subsys??_env")
os.system("rm -fr subsys?_env")

print ""
print ""
print " ---------------- script for making xcp working directory -------------"
print ""
print ""
print "files you need: "
print "    g.in_base   param.in   "


#--------- get information of this job ----------------
nspin         =  os.popen("grep \"nspin\" param.in | awk '{print $2}'").read()
total_mag     = float( os.popen("grep \"total_mag\" param.in | awk '{print $2}'").read() ) 
total_charge  = float( os.popen("grep \"total_charge\" param.in | awk '{print $2}'").read() ) 
nsubsys       =  os.popen("grep \"nsubsys\" param.in | awk '{print $2}'").read()
npspfile      =  os.popen("grep \"npsp\" param.in | awk '{print $2}'").read()
natom         = int( os.popen("grep \"natom\" param.in | awk '{print $2}'").read())
ntype         = int( os.popen("grep \"ntypat\" param.in | awk '{print $2}'").read() )
string_znucl  = os.popen("grep znucl param.in").read() 
part_method   =  os.popen("grep \"partition_method\" param.in | awk '{print $2}'").read()
tsmear        = float( os.popen("grep \"tsmear\" param.in | awk '{print $2}'").read() )
print "partition_method => ", part_method
part_method = part_method.strip()


# ------- get atomic numbers ----------
sp = string_znucl.split() 
z_list = []
for j in range(0,len(sp)): 
   ll = sp[j]
   # discard in the first keyword 
   if j==0: 
     continue 
   if ll[0] == '#': 
      break
   else:
     z_list.append(int(ll)) 

print "nspin: ",nspin
print "nsubsys: ",nsubsys
print "npspfile: ",npspfile 
print "z_list: ",z_list 


# ------------ get spin and charge -------------
spin = np.zeros(natom)
charge = np.zeros(natom)
if part_method != 'dm': 
   f = open("param.in",'r')
   lines = f.readlines()
   c = 0 
   for l in lines:
      l = l.strip()
      c = c+1
      if l == "subsys_charge_spin":
         print "\n\n ------ atom_ID ------- charge ------- spin --------\n"
         for q in range(natom):
            tmp = lines[c+q]
            sp = tmp.split()
            at = int(sp[0])
            charge[q] = float(sp[1])
            spin[q] = float(sp[2])
            print "atom: %d  charge: %f  spin: %f" %(at,charge[q],spin[q])
         break
   f.close()


# get acell  ---------------------------
########################################
f = open("param.in",'r')
lines = f.readlines()
c = 0
for l in lines:
   l = l.strip()
   if l == "cell_acell_angst":
      aa = lines[c+1]; aa= aa.strip()
      bb = lines[c+2]; bb= bb.strip()
      cc = lines[c+3]; cc= cc.strip()
      break
   c = c+1
f.close()      
print ""
print "cell information [read from param.in]"
print "  aa: ",aa
print "  bb: ",bb
print "  cc: ",cc
print ""
print ""


# ---------- get psps info ---------
print 'loading param.in ...'
f = open("param.in",'r')

pspf = []
counter  = 0
start_to_record = 0

for l in f.readlines():
  l=l.strip()
  if start_to_record == 1: 
     pspf.append( l ) 
     print pspf
     counter = counter + 1 
     if  counter == int(npspfile):
        break
  if l == "psps_file": 
     start_to_record = 1

f.close() 
print "psp files are ", pspf





# ------------- get cluster and env informations from param.in file -------------
atom_types = []
xangst = np.ones((natom,3))          # cooridnates of all atoms 
nbuffer = []                         # how many buffer atoms for each atom 
nbuffer_index = np.zeros((natom,20),dtype=np.int)   # what are the index of each buffer atom 
cluster_info_flag = -1 
counter = 0
f = open("param.in",'r')


for l in f.readlines(): 
    l = l.strip() 

    if l == "atom_info": 
        cluster_info_flag = 1
        continue

    #  read cluster info 
    #
    if cluster_info_flag == 1: 

        print "loading atom_info for",counter
        lsp = l.split()
        atom_types.append( int( lsp[1]) )
        xangst[counter][0] = float(lsp[2])
        xangst[counter][1] = float(lsp[3])
        xangst[counter][2] = float(lsp[4])

        # load the list of buffer atoms 
        if int(lsp[5])==999:
           nbuffer.append( natom-1 )
           pp = 0;
           for q in range(0,natom): 
              if q != counter:
                 nbuffer_index[counter][pp] = int( q+1 )
                 pp = pp+1
        else:        
           nbuffer.append(int(lsp[5]))
           # loop over buffer atom in that line 
           for q in range(0,nbuffer[counter]): 
              nbuffer_index[counter][q] = int(lsp[q+6] )

        counter = counter + 1
    
    # after reading all atom in the atom_info block, exit
    if counter >= natom:
       break
       
f.close()
print "atom types: \n", atom_types 
print "xangst: \n",xangst
print "nbuffer: \n",nbuffer 
print "nbuffer_index: \n",nbuffer_index
print "\n\n"


#
# 
# 
# ------- global system 
# 
# 
# 
print "\nmaking py_g.in_global file..... \n"
ff_gb = open("tmp_input_gb",'w')  # for global system
ff_gb.write("")
ff_gb.write("")
ff_gb.write("#--------- automatically generated by prepare_xcp.py --------- \n")
ff_gb.write("")
ff_gb.write("")
ff_gb.write("tsmear %f  eV\n"%(tsmear))
ff_gb.write("\nacell 1.0 1.0 1.0  angstrom\n")
ff_gb.write("rprim \n")
ff_gb.write("%s \n"%(aa))
ff_gb.write("%s \n"%(bb))
ff_gb.write("%s \n"%(cc))
ff_gb.write("\n \n")
ff_gb.write("ntypat  %i \n"%(ntype))
ff_gb.write(string_znucl)
string = "natom %i"%(int(nsubsys))
ff_gb.write(string)
ff_gb.write("\nxangst \n")

# get types and cooridnates of atoms in the total system
ty_index_global = []
for i in range(1,int(nsubsys)+1): 
   string = "   %f %f %f \n"%(xangst[i-1][0], xangst[i-1][1], xangst[i-1][2]) 
   ff_gb.write(string) 
   for p in range(0,ntype): 
      if atom_types[i-1]  == z_list[p]: 
         ty_index_global.append(p+1) # find the index of the atom type for this atom 
   
if int(nspin)==2:
   ff_gb.write("\n\nspinmagntarget %f "%(total_mag))
ff_gb.write("\n\ncharge         %f "%(total_charge))

# typat block 
ff_gb.write("\n\ntypat \n")
for u in range(0,len(ty_index_global)): 
   string = "% i"%(ty_index_global[u]) 
   ff_gb.write(string) 
ff_gb.write("\n\n")   

ff_gb.close()
cmd = "cp g.in_base py_g.in_global"; os.system(cmd) 
cmd = "cat tmp_input_gb  >> py_g.in_global"; os.system(cmd) 


#
# subsystem 
#
# --------------- make input files g.in for global system and subsystems -------------
print "make input files g.in files ............ \n"
for i in range(1,int(nsubsys)+1): 

   print '---- atom: ',i,' --------'

   ff = open("tmp_input",'w')
   if part_method != "dm": 
       ff_env = open("tmp_input_env",'w')

   ff.write("\n\n\n#--------- automatically generated by prepare_xcp.py --------- \n\n\n")
   if part_method !="dm": 
     ff_env.write("\n\n\n#--------- automatically generated by prepare_xcp.py --------- \n\n\n")

   ff.write("tsmear %f  eV\n"%(tsmear))
   ff_env.write("tsmear %f eV \n"%(tsmear))

   ff.write("\nacell 1.0 1.0 1.0  angstrom\n")
   ff.write("rprim \n")
   ff.write("%s \n"%(aa))
   ff.write("%s \n"%(bb))
   ff.write("%s \n"%(cc))
   ff.write("\n \n")
   if part_method != "dm": 
      ff_env.write("\nacell 1.0 1.0 1.0  angstrom\n")
      ff_env.write("rprim \n")
      ff_env.write("%s \n"%(aa))
      ff_env.write("%s \n"%(bb))
      ff_env.write("%s \n"%(cc))
      ff_env.write("\n \n")

   if int(nspin)==2: 
      ff.write("spinmagntarget  %f \n"%(spin[i-1]))
   ff.write("charge          %f \n"%(charge[i-1]))

   ff.write("ntypat  %i \n"%(ntype))
   if part_method!="dm": 
       if int(nspin)==2: 
         ff_env.write("spinmagntarget  %f \n"%(total_mag    - spin[i-1]))
       ff_env.write("charge          %f \n"%(total_charge - charge[i-1]))
       ff_env.write("ntypat  %i \n"%(ntype))

   ff.write(string_znucl)
   if  part_method!="dm": 
       ff_env.write(string_znucl)

   ty_index_cluster = []
   ty_index_env  = []
   natom_local = 0


   # ------ coordinates --------
   ff.write("\nxangst \n")
   ff_env.write("\nxangst \n")

   #
   # loop over all  atoms 
   #
   for m  in range(1,natom+1): 
      cluster_atom = -1
      if m==i: 
          ###### m is the central atom #######
          string = "   %f %f %f # center atom [index: %d]\n"%(xangst[i-1][0], xangst[i-1][1], xangst[i-1][2],m) 
          ff.write(string) 
          natom_local = natom_local + 1
          cluster_atom = 1
      else: 
          # check m is the buffer atom ?
          # loop over all buffer atoms of atom_i 
          for j in range(0,nbuffer[i-1]): 
             # qq is the index of the buffer atom 
             qq = nbuffer_index[i-1][j]-1
             if m-1==qq: 
               string = "   %f %f %f # buffer atom [index: %d]\n"%(xangst[qq][0],xangst[qq][1],xangst[qq][2],qq+1)  
               ff.write(string) 
               natom_local = natom_local + 1
               cluster_atom  = 1

      ####### m is the environmental atom #######
      if part_method!='dm' and cluster_atom==-1: 
          string = "   %f %f %f # env atoms \n"%(xangst[m-1][0], xangst[m-1][1], xangst[m-1][2]) 
          ff_env.write(string) 

      # ------- find the type index of this atom j --------
      ty = atom_types[m-1] 
      for p in range(0,ntype): 
         if ty == z_list[p]: 
             if cluster_atom == 1: 
                ty_index_cluster.append(p+1) 
             else: 
                ty_index_env.append(p+1) 
   
   ####################################################################
   ################ end of loop of all atoms in the cell ##############
   ####################################################################

   ff.write("\ntypat \n")
   for u in range(0,len(ty_index_cluster)): 
      string = "% i"%(ty_index_cluster[u]) 
      ff.write(string) 
   ff.write("\n\n")

   # -------------------------- 
   if part_method != 'dm': 
       ff_env.write("\ntypat \n")
       for u in range(0,len(ty_index_env)): 
          string = "% i"%(ty_index_env[u]) 
          ff_env.write(string) 
       ff_env.write("\n\n")

   #---------------  
   string = "natom %i"%(natom_local )
   ff.write(string)
   if part_method != 'dm': 
       string = "natom %i"%(natom - natom_local)
       ff_env.write(string)

   ff.close()
   if part_method != 'dm': 
       ff_env.close()

   cmd = "cp g.in_base  py_g.in_sub%i"%(i);  os.system(cmd) 
   cmd = "echo \# below are automatically generated by prepare_xcp.py ------- >> py_g.in_sub%i"%(i);  os.system(cmd) 
   cmd = "echo   >> py_g.in_sub%i"%(i);  os.system(cmd) 
   cmd = "echo   >> py_g.in_sub%i"%(i);  os.system(cmd) 
   cmd = "cat tmp_input  >> py_g.in_sub%i"%(i); os.system(cmd) 

   if part_method != 'dm': 
       cmd = "cp g.in_base         py_g.in_sub%i_env"%(i);  os.system(cmd) 
       cmd = "cat tmp_input_env  >> py_g.in_sub%i_env"%(i); os.system(cmd) 





# -------------- global system ----------------
print "\n removing global_system folder ........."
cmd = "rm -rf ./global_system "; os.system(cmd) 
print "making global_system folder .........\n"
cmd = "mkdir global_system";   os.system(cmd)
cmd = "cp g.files  global_system/"; os.system(cmd)
print cmd 
cmd = "cp py_g.in_global global_system/g.in"; os.system(cmd)
print cmd 




# ------------- preparing subsystem folders -----------
print "make subsystem folders ............ \n"

for i in range(1, int(nsubsys)+1): 

  cmd = "mkdir subsys%i"%(i);   os.system(cmd )
  if part_method !='dm':
        cmd = "mkdir subsys%i_env"%(i);   os.system(cmd )

  cmd = "cp py_g.in_sub%i subsys%i/g.in "%(i,i);  os.system(cmd) 
  print cmd
  cmd = "cp g.files  subsys%i/ "%(i); os.system(cmd) 
  if part_method !='dm': 
        cmd = "cp py_g.in_sub%i_env subsys%i_env/g.in "%(i,i);  os.system(cmd) 
        print cmd
        cmd = "cp g.files  subsys%i_env/ "%(i); os.system(cmd) 

  # ----------- copy psp files to folders ----------
  for q in range(0,  int(npspfile)):
     if part_method != 'dm': 
          cmd = "cp %s subsys%i_env/"%(pspf[q],i);  os.system(cmd)
     cmd = "cp %s subsys%i/    "%(pspf[q],i);  os.system(cmd)
     cmd = "cp %s global_system/ "%(pspf[q]);  os.system(cmd)
