# -*- coding: utf-8 -*-

"""
Temporal network computational experiments
"""
from scipy import sparse
from scipy import linalg
from sympy import Symbol
import networkx as nx
import numpy as np
from matplotlib.cm import ScalarMappable
import scipy.io
from pathlib import Path
import matplotlib.pyplot as plt

networkFileName = 'D:\\fauci-email-temporalgraph-tofrom'

savePath = "D:\\Matrices\\"

#Options
#binarise: Set all weights to 1.
#saveFile: Saves the matrices produced.
binarise = False
saveFile = False



#Read the data from the .json file.
with open(networkFileName + ".json", "r") as f:
  f.readline() # read the first '{'
  nverts = int(f.readline().split(':')[1].split(',')[0]) #read the number of vertices for the entire graph
  ngraphs = int(f.readline().split(':')[1].split(',')[0]) #read the number of graphs
  f.readline() # read the label "sequencedata"
  no_edges, edges, srcs, dsts, weights = [], [], [], [], []
  for _ in range(ngraphs):
    src, dst, weight = [], [],[]
    no_edges.append(int(f.readline().split(":")[1].split(',')[0]))
    f.readline() #edgedata_line 
    for i in range(no_edges[-1]): #iterative over the value of keyword "edges"
        einfo = f.readline().split(",") #split each edge by ','
        src.append(int(einfo[0]))
        dst.append(int(einfo[1]))
        weight.append(int(einfo[2]))
    srcs.append(src)
    dsts.append(dst)
    edges.append(list(zip(src, dst, weight))) 
    weights.append(weight)
    f.readline() # read end array
    f.readline() # read label array start
    f.readline() # read end graph start
 
    
  print('srcs, dsts, weights, edges, no_edges produced.')
  
  f.readline() #read the line '"dates": ['
  dates = []
  for _ in range(ngraphs):
    dates.append(f.readline().strip().strip(",").strip('"'))
  f.readline() # read dates array end
  print('dates produced')
  
  f.readline() # read labels array start
  labels = []
  for _ in range(nverts):
    labels.append(f.readline().strip().strip(","))
  f.readline() # read labels array end
  print('labels produced')
  
  
  f.readline()  # read orgs array start
  orgs = []
  for _ in range(nverts):
      orgs.append(int(f.readline().split(",")[0]))
  print('orgs produced')





#Next we loop over the number of graphs, producing at each time frame L,R, A, W, B and S.
Ls = []
Rs = []
As = []
Ss = []
Katzrhos = []
for t in range(ngraphs):
    L = np.zeros([no_edges[t], nverts])
    R = np.zeros([no_edges[t], nverts])
    S = []
    for e in range(no_edges[t]):
        i,j,w = edges[t][e]
        L[e][i] = 1
        R[e][j] = 1
        S.append(w)
    Ls.append(sparse.csr_matrix(L))
    Rs.append(sparse.csr_matrix(R))
    if binarise:
        As.append(sparse.csr_matrix(L.T@R))
    else:
        As.append(sparse.csr_matrix(L.T@np.diag(S)@R))
    Katzrhos.append(max(np.linalg.eig(L.T@np.diag(S)@R)[0], key=lambda x: abs(x)))
    Ss.append(sparse.csr_matrix(np.diag(S)))
    


#Now we make M, and also record the Perron eigenvalue of each diagonal block
rows = []
NBTrhos = []

for i in range(ngraphs):
    row = []
    for j in range(ngraphs):
        if (i <= j): #We don't set to zero
            W = Ss[i]@Rs[i]@(Ls[j].T)@Ss[j]
            W_time_reversed = Ss[j]@Rs[j]@(Ls[i].T)@Ss[i]
            
            
            backtracking_part = np.sqrt(W.multiply(W_time_reversed.T)) 
            
            
            B = W - backtracking_part
          
            B = np.sqrt(B)

            
            row.append(B)

            if (i == j):
                if not binarise:
                    NBTrhos.append(
                        max(sparse.linalg.eigs(B.toarray(),1)[0], key = lambda x: np.abs(x))
                            )
                else:
                    binarise_B = B.copy()
                    binarise_B.data[:] = np.ones([1, len(B.data)])
                    NBTrhos.append(
                        max(sparse.linalg.eigs(binarise_B.toarray(),1)[0], key=lambda x: np.abs(x))
                        )
                                    
        
        else: #We set to zero
            
            row.append(np.zeros([Ls[i].shape[0], Ls[j].shape[0]]))
    #print('Row {}/{}'.format(i+1, ngraphs))
    rows.append(sparse.hstack(row))
    
M = sparse.vstack(rows)
bigL = sparse.vstack(Ls)
bigR = sparse.vstack(Rs)
sqrtS =  np.sqrt(sparse.block_diag(Ss))
binarised_label = ''
if binarise:
    binarised_M = M.copy()
    binarised_M.data[:] = np.ones([1, len(M.data)])
    M = binarised_M
    del binarised_M
    binarised_label = ' binarised '

NBTrho = max(NBTrhos, key= lambda x: abs(x))
Katzrho = max(Katzrhos, key= lambda x: abs(x))

both_rhos = np.array([Katzrho,NBTrho])
    
if saveFile:
    scipy.io.savemat(savePath+ binarised_label + "Temporal_As"+ '.mat', mdict={'As'.format(i): As} )
    scipy.io.savemat(savePath+ binarised_label + "Temporal_sqrtS"+ '.mat', mdict={'sqrtS': sqrtS} )
    scipy.io.savemat(savePath+ binarised_label + "Temporal_bigL"+ '.mat', mdict={'bigL': bigL} )
    scipy.io.savemat(savePath+ binarised_label + "Temporal_bigR"+ '.mat', mdict={'bigR': bigR} )
    scipy.io.savemat(savePath + binarised_label +"Temporal_M"+ '.mat', mdict={'M': M} )
    scipy.io.savemat(savePath + binarised_label +"Temporal_rhos_KatzNBT"+ '.mat', mdict={'both_rhos': both_rhos} )
