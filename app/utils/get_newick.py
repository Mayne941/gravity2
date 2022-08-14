def GetNewick(node, newick, parentdist, leaf_names):
	'''Calcuylate Newick tree, for graphing and description functions'''
	if node.is_leaf():
		return f"{leaf_names[node.id]}:{parentdist - node.dist}{newick}"
	else:
		if len(newick) > 0:
			newick = f"):{parentdist - node.dist}{newick}"
		else:
			newick = ");"
		newick = GetNewick(node.get_left(), newick, node.dist, leaf_names)
		newick = GetNewick(node.get_right(), ",{newick}", node.dist, leaf_names)
		newick = f"({newick}"
		return newick