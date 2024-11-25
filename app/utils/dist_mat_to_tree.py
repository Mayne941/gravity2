import numpy as np
from scipy.cluster.hierarchy import linkage, to_tree
from .get_newick import GetNewick

def DistMat2Tree(DistMat, LeafList, Dendrogram_LinkageMethod):
	'''Convert (similarity) distance matrix to tree, for graphing functions'''
	DistMat[DistMat<0]	= 0
	DistList			= DistMat[np.triu_indices_from(DistMat, k = 1)]
	DO_LOGSCALE = False # RM < TODO Parameterise
	if DO_LOGSCALE:
		DistList = 10**DistList
		# DistList = (DistList-np.min(DistList))/(np.max(DistList)-np.min(DistList)) # Standard scaler
	linkageMat			= linkage(DistList, method = Dendrogram_LinkageMethod)
	TreeNewick			= to_tree(Z	 = linkageMat,
								  rd = False,
					)
	TreeNewick			= GetNewick(node		= TreeNewick,
									newick		= "",
									parentdist	= TreeNewick.dist,
									leaf_names	= LeafList,
					)
	return TreeNewick
