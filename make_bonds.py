#!/usr/bin/python

import numpy as np
import sys



##################### Dijkstra's Algorithm #################
# http://www.gilles-bertrand.com/2014/03/dijkstra-algorithm-python-example-source-code-shortest-path.html

def dijkstra(graph,src,dest,visited=[],distances={},predecessors={},print_out=None):
    """ calculates a shortest path tree routed in src
    """    
    # a few sanity checks
    if src not in graph:
        raise TypeError('The root of the shortest path tree cannot be found')
    if dest not in graph:
        raise TypeError('The target of the shortest path cannot be found')    
    # ending condition
    if src == dest:
        # We build the shortest path and display it
        path=[]
        pred=dest
        while pred != None:
            path.append(pred)
            pred=predecessors.get(pred,None)
        if print_out>0:
           print('shortest path: '+str(path)+" cost="+str(distances[dest])) 
        return str(distances[dest])
    else :     
        # if it is the initial  run, initializes the cost
        if not visited: 
            distances[src]=0
        # visit the neighbors
        for neighbor in graph[src] :
            if neighbor not in visited:
                new_distance = distances[src] + graph[src][neighbor]
                if new_distance < distances.get(neighbor,float('inf')):
                    distances[neighbor] = new_distance
                    predecessors[neighbor] = src
        # mark as visited
        visited.append(src)
        # now that all neighbors have been visited: recurse                         
        # select the non visited node with lowest distance 'x'
        # run Dijskstra with src='x'
        unvisited={}
        for k in graph:
            if k not in visited:
                unvisited[k] = distances.get(k,float('inf'))        
        x=min(unvisited, key=unvisited.get)
        return dijkstra(graph,x,dest,visited,distances,predecessors,print_out)


####################################
# example for using dijkstra()
#
#graph = {'s': {'a': 2, 'b': 1},
#            'a': {'s': 3, 'b': 4, 'c':8},
#            'b': {'s': 4, 'a': 2, 'd': 2},
#            'c': {'a': 2, 'd': 7, 't': 4},
#            'd': {'b': 1, 'c': 11, 't': 5},
#            't': {'c': 3, 'd': 5}}
#dijkstra(graph,'t','s')


################# load XYZ file ############

nb = int(raw_input("give Nb ==>  "))

f = open("coord.xyz", "r")

at_index = np.array([]);
at_type  = np.array([]);
at_x = np.array([]);
at_y = np.array([]);
at_z = np.array([]);

natom = 0

for line in f:
    natom = natom + 1;
    l = line.split()
    at_index = np.append(at_index,int(l[0]))
    at_type  = np.append(at_type,int(l[1]))
    at_x     = np.append(at_x,float(l[2]))
    at_y     = np.append(at_y,float(l[3]))
    at_z     = np.append(at_z,float(l[4]))

#print at_index
#print at_type 
#print at_x
#print at_y
#print at_z



################## make edges #################

nbond = 0
edges = {}


# threshold of the bonds 
bond_CH = 1.5
bond_OH = 1.5
bond_CC = 2.0
bond_HH = 1.3
bond_CO = 2.0
bond_OO = 2.1
bond_NO = 2.0
bond_CN = 2.0
bond_NH = 1.6
bond_NN = 2.0


for i in range(natom): 
   x1 = at_x[i]
   y1 = at_y[i]
   z1 = at_z[i]
   e1 = "%i"%(at_index[i])

   d = {};
   for j in range(natom): 
     if j==i: continue 

     x2 = at_x[j]
     y2 = at_y[j]
     z2 = at_z[j]

     dist = ((x1-x2)**2 + (y1-y2)**2 + (z1-z2)**2)**0.5

     ########### threshold on bond lengths ###########
     if (at_type[i]==1 and at_type[j]==1):
        threshold = bond_HH 
     elif (at_type[i]==6 and at_type[j]==6):
        threshold = bond_CC
     elif (at_type[i]==8 and at_type[j]==8):
        threshold = bond_OO
     elif (at_type[i]==6 and at_type[j]==8): 
        threshold = bond_CO
     elif (at_type[i]==8 and at_type[j]==6): 
        threshold = bond_CO
     elif (at_type[i]==6 and at_type[j]==1): 
        threshold = bond_CH
     elif (at_type[i]==1 and at_type[j]==6): 
        threshold = bond_CH
     elif (at_type[i]==8 and at_type[j]==1): 
        threshold = bond_OH
     elif (at_type[i]==1 and at_type[j]==8): 
        threshold = bond_OH
     else:
        print "bond not defined. stop"
        print at_type[i]
        print at_type[j]
        stop

     if (dist<threshold): 
        nbond = nbond+1
        e2 = "%i"%(at_index[j])
        d[e2] = 1
   
   edges[e1] = d


print "number of atoms: ",natom
print "number of bonds: ",nbond/2
print "edges: \n",edges
 
#print dijkstra(edges, '1', '25')
#stop

nlist = np.zeros([natom,natom])
nat = np.zeros([natom,1])

################ get the neighours ###########
for i in range(natom): 

   n = 0
   for j in range(natom): 

      if (j==i): continue 

      visited=[]
      distances={}
      predecessors={}
      print_out = -1
      node1 = str(int(at_index[i]))
      node2 = str(int(at_index[j]))
      cost = int(dijkstra(edges, node1, node2, visited, distances, predecessors, print_out))

      if cost<=nb: 
         nlist[i][n] = int(j+1)
         n = n + 1
         nat[i] = nat[i]+1


 
#print "\n\n\nnumber of atoms: ",natom
#print "number of bonds: ",nbond/2
#print "edges: \n",edges



############# print neighour list #########
print "\n\n ------------- neighbor list --------------\n"
print " # of atoms |   <atom list >"
for i in range(natom): 
   print "%4d        "%nat[i],
   for j in range(natom): 
      aa = nlist[i][j]
      if (aa!=0): 
         print "%3d"%int(aa),
   
   print ""

