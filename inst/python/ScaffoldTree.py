#------Libraries needed for Tree Topology Functions ------
import scipy as sp
import numpy as np
import collections
from scipy import stats
import sys
import csgraph_mod
from csgraph_mod import shortest_path_mod as spm
import time

#--------------Functions-------------------------------------
def SinglePath(DijkstraPredecesors, i, j):
    path = []
    k=j
    while (i != k) & (k >= 0):
        path.append(k)
        k = DijkstraPredecesors[i,k]
    path.append(i)
    return (path)

def calculate_cell_paths(Distancias, DijkstraPredecesors):
    Paths = collections.defaultdict(lambda: collections.defaultdict(list))
    for i in range(0, Distancias.shape[0]):
        for j in range(0, Distancias.shape[1]):
            path = []
            k=j
            while (i != k) & (k >= 0):
                path.append(k)
                k = DijkstraPredecesors[i,k]
            path.append(i)
            Paths[i][j]=path
    return(Paths)

def calc_extreme_endpoints(CellDistances):
    EndPoints=[]
    Endpoint1=np.where(CellDistances==np.max(CellDistances))[0][0]
    Endpoint2=np.where(CellDistances==np.max(CellDistances))[1][0]
    EndPoints.append(Endpoint1)
    EndPoints.append(Endpoint2)
    return EndPoints
    
def calculate_secondary_endpoints(Endpoints, DijkstraMatrix, NumberNodes, NBranches=-1, graph='no'):
    R_epsilon_todos={}
    TryEndpoint=True
    NewBranch=True
    NBranches=int(NBranches)
    nodes=list(range(DijkstraMatrix.shape[0]))
    BranchMinLength=np.sqrt(len(nodes)/2)
    while TryEndpoint:
        #Add a new endpoint
        ScoreEndpoints=np.zeros(len(nodes))
        ScoreEndpointsNodes=np.zeros(len(nodes))
        for n in nodes:
            if n in EndPoints:
                continue
            else:
                CombinatoryEndpointIncreament=[]
                CombinatoryEndpointIncreamentNodes=[]
                for k in range(0, len(EndPoints)):
                    for l in range(0, len(EndPoints)):
                        if k < l:
                            AddedDistance= DijkstraMatrix[n][EndPoints[k]] + DijkstraMatrix[n][EndPoints[l]] - DijkstraMatrix[EndPoints[k]][EndPoints[l]]
                            AddedDistanceNodes= NumberNodes[n][EndPoints[k]] + NumberNodes[n][EndPoints[l]] - NumberNodes[EndPoints[k]][EndPoints[l]]
                            CombinatoryEndpointIncreament.append(AddedDistance)
                            CombinatoryEndpointIncreamentNodes.append(AddedDistanceNodes)
                ScoreEndpoints[n]=0.5*min(CombinatoryEndpointIncreament)
                ScoreEndpointsNodes[n]=0.5*min(CombinatoryEndpointIncreamentNodes)

        #we get the nodes with the largest number of cells on path
        nodes_maxnodes=np.where(ScoreEndpointsNodes==np.max(ScoreEndpointsNodes))[0]

        #from the array of nodes with the largest number of cells on path, we take the one with the max shortest path
        new_endpoint=nodes_maxnodes[np.where(ScoreEndpoints[nodes_maxnodes]==np.max(ScoreEndpoints[nodes_maxnodes]))[0][0]]
        
        R_epsilon={}
        R_epsilon_val=[]

        for i in nodes:
            if i not in EndPoints and ScoreEndpointsNodes[i] > 0 :
                R_epsilon[i]=ScoreEndpointsNodes[i]
                R_epsilon_val.append(ScoreEndpointsNodes[i])
        
        if len (R_epsilon_val) >1:
            mean_epsilon=np.mean(R_epsilon_val)
            std_epsilon=np.std(R_epsilon_val)
            var_epsilon=np.var(R_epsilon_val)

        if NBranches > 0 :
        	if (len(EndPoints)) < NBranches:
        		EndPoints.append(new_endpoint)
        		if graph == "yes":
        			print_hist_newpoint(R_epsilon_val, R_epsilon, new_endpoint)
        	else:
        		TryEndpoint=False
        else:
        	if ScoreEndpointsNodes[new_endpoint]<BranchMinLength:
		        NewBranch=False

        	if R_epsilon[new_endpoint] >= BranchMinLength:
		        EndPoints.append(new_endpoint)
		        if graph == "yes":
		            print_hist_newpoint(R_epsilon_val, R_epsilon, new_endpoint)

        	else:
		        TryEndpoint=False

    return (R_epsilon_todos)

def calculate_branchpoints(EndPoints, DijkstraMatrix, DijkstraPredecesors, NumberNodes):
    NJDistances = {}
    VBranchpoints=EndPoints
    #selecting branching points
    Branchpoints=[]
    TreeConnectivity=[]

    while len(VBranchpoints) >2:
        NJDistances = {}
        for branchpoint1 in VBranchpoints:
            for branchpoint2 in VBranchpoints:
                if branchpoint1 < branchpoint2:
                    dbranchpoint1=0
                    dbranchpoint2=0

                    for m in VBranchpoints:
                        dbranchpoint1+=NumberNodes[branchpoint1][m]
                        dbranchpoint2+=NumberNodes[branchpoint2][m]
                    NJDistances[(branchpoint1,branchpoint2)]=(len(VBranchpoints)-2)*NumberNodes[branchpoint1][branchpoint2]-dbranchpoint1-dbranchpoint2
        
        minbranchpoint1=min(zip(NJDistances.values(), NJDistances.keys()))[1][0]
        minbranchpoint2=min(zip(NJDistances.values(), NJDistances.keys()))[1][1]
        
        mDistances={}
        for m in SinglePath(DijkstraPredecesors, minbranchpoint1, minbranchpoint2):
            if(m not in EndPoints):
                otherbranchpoints=list(set(VBranchpoints)-set([minbranchpoint1, minbranchpoint2]))
                mDistanceothers=0
                for o in otherbranchpoints:
                    mDistanceothers+=NumberNodes[o][m]
                mDistances[m]=NumberNodes[m][minbranchpoint1]+NumberNodes[m][minbranchpoint2]+(1/(len(VBranchpoints)-2))*mDistanceothers

        branchpoint=min(zip(mDistances.values(), mDistances.keys()))[1]
        VBranchpoints=list(set(VBranchpoints)-set([minbranchpoint1, minbranchpoint2]))

        TreeConnectivity.append([minbranchpoint1, branchpoint])
        TreeConnectivity.append([minbranchpoint2, branchpoint])

        VBranchpoints.append(branchpoint)
        Branchpoints.append(branchpoint)

		#Here we add as a new connection the one corresponding to the root and the last detected branchpoint.
        if len(VBranchpoints)==2:
        	TreeConnectivity.append([VBranchpoints[0], VBranchpoints[1]])
        	TreeConnectivity=np.array(TreeConnectivity)
    
    return Branchpoints, TreeConnectivity

def euclidean(v1, v2):
    distance=sp.spatial.distance.euclidean(v1, v2)
    return distance

def calc_distances(X):
    distances = np.zeros((X.shape[0], X.shape[0]))
    
    for i in range(X.shape[0]):
        for j in range(X.shape[0]):
            if (i < j):
                d = euclidean(X[i], X[j])
                distances[i][j] = d
                distances[j][i] = d
    return(distances)

def LengthPath(DijkstraPredecesors, i, j):
    length=0
    k=j
    while (i != k) & (k >= 0):
        length+=1
        k = DijkstraPredecesors[i,k]
    return (length)

def calculate_bi_longerpath(Distancias, DijkstraPredecesors):
    NumberNodes=np.zeros(Distancias.shape)
    first=-1
    second=-1
    ncells=0
    dist=-1
    
    for i in range(0, Distancias.shape[0]):
        if(i%100==0):
            print("cell", i, "out of", Distances.shape[0], "analized...")
        for j in range(0, Distancias.shape[1]):
            if(i<j):
                path=LengthPath(DijkstraPredecesors, i, j)
                if(path>ncells and Distancias[i][j]>dist):
                    first=i
                    second=j
                    dist=Distancias[i][j]
                    ncells=path
                NumberNodes[i][j]=path
                NumberNodes[j][i]=path
    return(first, second, NumberNodes)

print ("Loaded Tree Topology Library...")

#----------------------------------------------------------------
#-----------------Main Program-----------------------------------
#----------------------------------------------------------------

import scipy.linalg as linalg
import os, math
import random
import sys
from scipy.sparse import csgraph
import pandas as pd
import argparse

# Arguments Parser and help
parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('Filename', help='Matrix with "cells" x "genes" dimensions. Fields shout be delimited by tabs')
parser.add_argument('-NBranches', type=int, default=(-1), help='Number of desired branches for the resulting tree. By default the value is set to -1, meaning that all branches longer than sqrt(N/2) will be considered')
parser.add_argument('-showplot', type=str, default="no", help='Weather or not 2D/3D plots should be shown. Options are "yes" and "no". By default it is not activated')
args = parser.parse_args()
DMCoordinates=args.Filename
NBranches=args.NBranches
plot=args.showplot

#---Read Manifold Coordinates (DMs, t-SNE, etc)
Coordinates = pd.read_csv(DMCoordinates, sep = "\t", header=None)
Coordinates=np.array(Coordinates)
DataDimensions=Coordinates.shape
print("Data Dimensions", DataDimensions)

start_distances = time.time()
print("Calculating distances...")
Distances=calc_distances(Coordinates)
InputMatrix=np.power(Distances, 2)

end_distances=time.time()
print("Calculating Scaffold Tree...")

#Calculate Shortest Distances and Predecesors
start_dijkstra = time.time()
DijkstraMatrix, DijkstraSteps, DijkstraPredecesors = spm.shortest_path(InputMatrix, method="D", return_predecessors=True, directed=False)
DijkstraSteps=DijkstraSteps-1

#Calculate Endpoints in the tree
EndPoints=calc_extreme_endpoints(DijkstraSteps)
Epsilons=calculate_secondary_endpoints(EndPoints, DijkstraMatrix, DijkstraSteps, NBranches)

if len(EndPoints) <= 2: 
	print ("Only 2 Endpoints detected. There is no tree to reconstruct.")
	print("Endpoints:", len(EndPoints), *EndPoints)

	TreeTopologyDat=DMCoordinates+"_TreeTopology.dat"
	with open(TreeTopologyDat, 'w') as out:
		print("Endpoints:", len(EndPoints), *EndPoints)
		EndPointsLine="Endpoints:"+"\t"+str(len(EndPoints))+"\t"+ ' '.join(map(str, EndPoints))+"\n"
		out.write(EndPointsLine)
	out.close()

	DijkstraStepsDat=DMCoordinates+"_DijkstraSteps.dat"
	DijkstraDistancesDat=DMCoordinates+"_DijkstraDistances.dat"
	DijkstraPredecesorsDat=DMCoordinates+"_DijkstraPredecesors.dat"

	np.savetxt(DijkstraDistancesDat, DijkstraMatrix)
	np.savetxt(DijkstraStepsDat, DijkstraSteps)
	np.savetxt(DijkstraPredecesorsDat, DijkstraPredecesors)

else:
	#Calculate branchpoints and tree topology
	branching, TreeConnectivity =calculate_branchpoints(EndPoints, DijkstraMatrix, DijkstraPredecesors, DijkstraSteps)
	print("Endpoints:", len(EndPoints), *EndPoints)

	TreeTopologyDat=DMCoordinates+"_TreeTopology.dat"
	with open(TreeTopologyDat, 'w') as out:
		print("Endpoints:", len(EndPoints), *EndPoints)
		EndPointsLine="Endpoints:"+"\t"+str(len(EndPoints))+"\t"+ ' '.join(map(str, EndPoints))+"\n"
		out.write(EndPointsLine)

		branching=set(branching)
		print("Branchpoints:", len(branching), *branching)
		BranchpointsLine="Branchpoints:"+"\t"+str(len(branching))+"\t"+' '.join(map(str, branching))+"\n"
		out.write(BranchpointsLine)

		#Shows tree connectivity:
		for i in range(TreeConnectivity.shape[0]):
			# In case of trifurcations or higher order connections it avoids printing branches in between the same node
			if TreeConnectivity[i][0] != TreeConnectivity[i][1]:
				print("Tree_Branch", TreeConnectivity[i][0], TreeConnectivity[i][1])
				TopologyLine="Tree_Branch"+"\t"+str(TreeConnectivity[i][0])+"\t"+str(TreeConnectivity[i][1])+"\n"
				out.write(TopologyLine)
	out.close()

	DijkstraStepsDat=DMCoordinates+"_DijkstraSteps.dat"
	DijkstraDistancesDat=DMCoordinates+"_DijkstraDistances.dat"
	DijkstraPredecesorsDat=DMCoordinates+"_DijkstraPredecesors.dat"

	np.savetxt(DijkstraDistancesDat, DijkstraMatrix)
	np.savetxt(DijkstraStepsDat, DijkstraSteps)
	np.savetxt(DijkstraPredecesorsDat, DijkstraPredecesors)

end_dijkstra = time.time()
print("Finished Scaffold Tree...", end_dijkstra-start_dijkstra, "seconds")
