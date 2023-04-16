'''
Shannon entropy is a way to calculate the probability of finding two similar objects together.

$$ S = - \frac{1}{N} \sum_i^N \sum_j^M \rho_j(i) \log \rho_j(i) $$
 where $\rho_j(i)$ is the number of particles of type $j$ inside $i$ divided by the total
 number of particles in $i$. The entropy of mixing of the whole
 system is then estimated as the average of the entropies over all
 subregions ($N$). 

Following code is a collection functions to calculate the Shannon entropy
between a two component system, whose residues are specified as "VT30" and
""P8EL" in GROMACS structure files (.gro). The code requires the .gro files
that are the snapshots of the system at regular time-steps where only the O
atom co-ordinates of the system are extracted.

For details please read: https://doi.org/10.1016/j.molliq.2023.121803

'''
import numpy as np
import pandas as pd
import os
import collections


def compute_entropy(fin="P8EL", # input file
                   n=3, # number of subregions
                   ):
   
    #--- check if in the box
    def check_if_in_range(point,pair,C):
        if point >= C[pair[0]] and point < C[pair[1]] :
            return True
        else:
            return False
    #-- get the probability
    def rho(box_contents,typ):
        return box_contents[typ]/sum(box_contents)
   
    #-- parse
    delta=0.02 # use to slightly increase the unit cell to accomodate some stray atoms
    f=open(fin,'r')
    data=f.readlines()
    unit_cell=[float(i)+delta for i in data[-1].split()]
    n_atms=int(data[1]) #number of atoms
    cords=data[2:-1] # cords
    if len(cords) != n_atms : #sanity check
        print (f"ERROR! There are {n_atms} atoms, but {len(cords)} cordinates.")
    f.close()
   
    #generate information of atoms -- type and positions
    atom_dict={}
    for c in cords:
        if 'VT30' in c.split()[0] :
            typ='VT30'
        else:
            typ='P8EL'
        pos=[float(j) for j in c.split()[-3:]]
        pos.append(typ)
        atom_dict[c.split()[0]]=pos
       
    #-- pairs for vertices
    X=[i*(unit_cell[0]/n) for i in range(n+1)]
    Y=[i*(unit_cell[1]/n) for i in range(n+1)]
    Z=[i*(unit_cell[2]/n) for i in range(n+1)]
   
    pairs=[(j,j+1) for j in range(n)]
   
    #-- generate boxes
   
    boxes=[]
    for i in pairs : # for X
        for j in pairs : # for Y'
            for k in pairs : # for Z
                boxes.append((i,j,k))
    #print (f"Total of {len(boxes)} regions")
   
   
    box_dict={} # what box contains what
    for key in atom_dict:
        check_point=atom_dict[key][:3]
        found=False

        for box in range(len(boxes)):
            if  check_if_in_range(check_point[0],boxes[box][0],X) \
            and check_if_in_range(check_point[1],boxes[box][1],Y) \
            and check_if_in_range(check_point[2],boxes[box][2],Z) :
                #print (f"Value {check_point[0]}, box number {box,} , X={X[boxes[box][0][0]]}, {X[boxes[box][0][1]]}" )
                #print (f"I'm box number {box} and I HAVE your point! X: {X[boxes[box][0][0]]:.1f} -> {X[boxes[box][0][1]]:.1f}, Y: {Y[boxes[box][1][0]]:.1f} -> {Y[boxes[box][1][1]]:.1f}, Z={Z[boxes[box][2][0]]:.1f} -> {Z[boxes[box][2][1]]:.1f}" )
                found=True
                if box in box_dict.keys():
                    box_contains=box_dict[box]
                    box_contains.append(atom_dict[key][-1])
                    box_dict[box]=box_contains
                else:
                    box_dict[box]=[atom_dict[key][-1]]

        if not found:
            print (f"{check_point} not found {unit_cell} ")
           
    box_freq_dict={}
    for box in range(len(boxes)):
        box_freq_dict[box]=[0,0] # (vt30,p8el)
        if box in box_dict.keys():
            frequency = dict(collections.Counter(box_dict[box]))
            #print (frequency)
            if 'VT30' in frequency.keys():
                c=box_freq_dict[box]
                c[0]=frequency['VT30']
                box_freq_dict[box]=c
            if 'P8EL' in frequency.keys():
                c=box_freq_dict[box]
                c[1]=frequency['P8EL']
                box_freq_dict[box]=c
    #sanity chk
    tolt_atoms=0
    for box in range(len(boxes)):
        tolt_atoms+=sum(box_freq_dict[box])
    if tolt_atoms != n_atms:
        print ("Err! Sanity check failed. Some atoms are not in any box?")
   
    # main calculation of entropy--
    rho_sum=0
    for box in range(len(boxes)):
        for j in range(2): #two types 0 for VT30, 1 for P8EL
            #print (box,j,box_freq_dict[box],rho(box_freq_dict[box],j))
            if sum(box_freq_dict[box]) != 0 :
                rho_ij=rho(box_freq_dict[box],j)
                if rho_ij != 0 : #dont take log(0)
                    rho_sum+=rho_ij*np.log(rho_ij)


    return -rho_sum/len(boxes) #entropy



def entropy4trajectory(dirname):
    # filenames are assummed to be O_{i}.gro .. 2000 files
    traj_ent=[]
    for i in range(1500,2001):
        fin=f"{dirname}/O_{i}.gro"
        if not os.path.isfile(fin): 
            print(f"File {dirname}/{fin} does not exist.")
        else:
            traj_ent.append(compute_entropy(fin,n=4))
    return traj_ent



def entropy4trajectorymono(dirname):
    traj_ent=[]
    for i in range(1500,2001):
        fin=f"{dirname}/O_mono_{i}.gro"
        if not os.path.isfile(fin): 
            print(f"File {dirname}/{fin} does not exist.")
        else:
            traj_ent.append(compute_entropy(fin,n=4))
    return traj_ent

def entropy4trajectoryaway(dirname):
    traj_ent=[]
    for i in range(1500,2001):
        fin=f"{dirname}/O_away_{i}.gro"
        if not os.path.isfile(fin): 
            print(f"File {dirname}/{fin} does not exist.")
        else:
            traj_ent.append(compute_entropy(fin,n=4))
    return traj_ent


