# -*- coding: utf-8 -*-
"""
Created on Mon Oct 17 14:39:29 2022

@author: rzmwo
"""

# -*- coding: utf-8 -*-
"""
Read the JSON file.
"""
import scipy as sp
import networkx as nx
import numpy as np
import scipy.io



savePath = 'D:\\Fauci paper redo\\Static Networks\\Hypergraph-projection-cc\\Matrices\\'
networkFileName = 'D:\\Fauci paper redo\\Static Networks\\Hypergraph-projection-cc\\json data\\fauci-email-graph-hypergraph-projection-cc'


saveMatrices = True

binarised = True;



"""
Read the networkFileName and extract the information
"""

"""
The following code produces:
    src, dst, weights

"""
with open(networkFileName + ".json", "r") as f:
  f.readline() # read the first '{'
  nverts = int(f.readline().split(':')[1].split(',')[0])
  nedges = int(f.readline().split(':')[1].split(',')[0])
  f.readline() # read "edgedata"
  src, dst, weights = [],[],[]
  for _ in range(nedges):
    einfo = f.readline().split(",")
    src.append(int(einfo[0]))
    dst.append(int(einfo[1]))
    weights.append(int(einfo[2]))
  f.readline() # read end array
  f.readline() # read label array start
  labels = []
  for _ in range(nverts):
    label_uncapitalized = f.readline().strip().strip(",").strip('"').split(',')
    label_capitalized= ', '.join([(k.strip(' ')).capitalize() for k in label_uncapitalized])
    labels.append(label_capitalized)
  f.readline() # read label array end
  f.readline() # read org array start
  orgs = []
  for _ in range(nverts):
    orgs.append(int(f.readline().strip().strip(",")))
    
    

    
#Then we produce the DiGraph with the edges read above.
graph = nx.DiGraph()
graph.add_weighted_edges_from(list(zip(src,dst,weights)))





""" Set name and organisation attributes """
nx.set_node_attributes(graph, {i: labels[i] for i in range(nverts)}, 'Name')
nx.set_node_attributes(graph, {i: orgs[i] for i in range(nverts)}, 'Organisation')





"""Produce the weighted adjancency matrix"""
A = nx.adjacency_matrix(graph, nodelist = list(range(nverts)))
m,n = A.shape




"""
NBT matrix M, L, R with edge order as in JSON file
"""
M = np.zeros([nedges, nedges])
SourceL = np.zeros([nedges, nverts])
TargetR = np.zeros([nedges, nverts])
sqrtS = scipy.sparse.csr_matrix(np.diag(weights)**0.5)

edges = list(zip(src, dst))
for i in range(len(edges)):
    u,v = edges[i]
    TargetR[i][v] = 1
    SourceL[i][u] = 1
    for j in range(len(edges)):
        r,s = edges[j]
        if (u != s) and (r == v):
            M[i][j] = weights[i]**0.5 *weights[j]**0.5
        else:
            continue

M = sp.sparse.csr_matrix(M)

binarised_label = ''
# If true then all weights are set to 1.
if binarised:
    binarised_label = 'binarised '
    Mcopy = M.copy()
    M.data[:] =np.ones([1,len(Mcopy.data[:])])
    sqrtS =np.eye(sqrtS.shape[0])
    A.data[:] = np.ones([1,len(A.data[:])])
    
"""
Rho the upper-limit of NBT Katz
"""

DNBT,VNBT = sp.sparse.linalg.eigs(M)
NBTrho =  max(DNBT, key= lambda x: abs(x))

D,V = np.linalg.eig(A.todense())
rho = max(D, key= lambda x: abs(x))

both_rhos = np.array([rho,NBTrho])
    

labels_numpy_obj = np.zeros((len(labels),), dtype=np.object)
labels_numpy_obj[:] = labels









"""
Now save the Katz/NBT vectors at 0.95/0.5 * 1/rho to txt files
"""
if binarised:
    isBinarised = 'Binarised'
else:
    isBinarised = ''

if saveMatrices:
    scipy.io.savemat(savePath+ binarised_label + "A"+ '.mat', mdict={'A': A} )
    scipy.io.savemat(savePath+ binarised_label + "sqrtS"+ '.mat', mdict={'sqrtS': sqrtS} )
    scipy.io.savemat(savePath+ binarised_label + "SourceL"+ '.mat', mdict={'SourceL': SourceL} )
    scipy.io.savemat(savePath+ binarised_label + "TargetR"+ '.mat', mdict={'TargetR': TargetR} )
    scipy.io.savemat(savePath + binarised_label +"M"+ '.mat', mdict={'M': M} )
    scipy.io.savemat(savePath + binarised_label +"rhos_KatzNBT"+ '.mat', mdict={'both_rhos': both_rhos} )
    scipy.io.savemat(savePath + binarised_label +"labels"+ '.mat', mdict={'labels': labels_numpy_obj} )