import numpy as np
from scipy.cluster.hierarchy import linkage, to_tree, cophenet
from scipy.spatial.distance import pdist

from .get_newick import GetNewick

def NormalizeData(data):
    return (data - np.min(data)) / (np.max(data) - np.min(data))

def DistMat2Tree(DistMat, LeafList, Dendrogram_LinkageMethod, do_logscale):
	'''Convert (similarity) distance matrix to tree, for graphing functions'''
	# DistMat[DistMat<0]	= 0 # SHOULD WE NORMALISE TO 0-1 WITHOUT THIS?
	# DistList			= DistMat[np.triu_indices_from(DistMat, k = 1)] # make squareform

	'''TEST'''
	DistMat_Norm = NormalizeData(DistMat)
	DistList	 = DistMat_Norm[np.triu_indices_from(DistMat_Norm, k = 1)] 
	''''''

	if do_logscale:
		DistList = 10**DistList

	# TODO < WIP: automatically determine best metric for dendrogram
	# options = {
	# 	"euclidean": [-1, None],
	# 	"cosine": [-1, None],
	# 	"cityblock": [-1, None],
	# 	"braycurtis": [-1, None],
	# 	"chebyshev": [-1, None],
	# }
	# options_idx_map = {key: i for key, i in enumerate(options)}

	# for option in options.keys():
	# 	linkageMat		= linkage(DistList, method = Dendrogram_LinkageMethod, metric=option, optimal_ordering=True) 
	# 	c_score, coph_dists = cophenet(linkageMat, pdist(DistMat_Norm, metric=option))
	# 	options[option][0] = c_score
	# 	options[option][1] = linkageMat

	# best_metric = options_idx_map[np.argmax([i[0] for i in options.values()])]
	# linkageMat = options[best_metric][1]

	linkageMat		= linkage(DistList, method = Dendrogram_LinkageMethod, metric='chebyshev', optimal_ordering=True) 
	TreeNewick			= to_tree(Z	 = linkageMat,
								  rd = False,
					)
	TreeNewick			= GetNewick(node		= TreeNewick,
									newick		= "",
									parentdist	= TreeNewick.dist,
									leaf_names	= LeafList,
					)
	return TreeNewick
