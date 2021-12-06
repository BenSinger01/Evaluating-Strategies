import numpy as np
from scipy import sparse
import networkx as nx

def random_walk_percolation_centrality(A,source):
	n = A.shape[0]
	n_source = len(source)
	PC_weighted = np.zeros(n)

	L_tilde = sparse.csgraph.laplacian(A)[1:n,1:n]
	L_tilde_inverse = np.linalg.inv(L_tilde)
	C = np.zeros((n,n))
	C[1:n,1:n] = L_tilde_inverse

	Ax = nx.from_numpy_matrix(A,create_using=nx.DiGraph)
	edges = list(Ax.edges)
	edge_by_vertex = [[i for i in range(len(edges)) if j in edges[i]] for j in range(n)]
	B = nx.linalg.graphmatrix.incidence_matrix(Ax,weight="weight",oriented=True)

	F = np.transpose(B)*C

	for v in np.delete(np.arange(n),source):
		for r in np.delete(np.arange(n),np.append(source,v)):
			for s in source:
				for edge in edge_by_vertex[v]:
					PC_weighted[v] += 0.5 * np.abs(F[edge,s]-F[edge,r]) / n_source
	return(PC_weighted)