#!/usr/bin/env python
import argparse
import os
import sys
import re
import colorsys
from ete3 import Tree, PhyloTree, faces, TreeStyle, RectFace, TextFace, add_face_to_node, AttrFace

# Adapted from http://etetoolkit.org/docs/latest/tutorial/tutorial_drawing.html#customizing-the-aspect-of-trees
# 	Image customization is performed through 4 main elements:
# 		1. Tree style:  controls the general aspect of the tree image
#		2. Node style:  controls various aspects of each single node
#		3. Node face:   graphical information that can be linked to nodes
#		4. Layout func: act as pre-drawing hooking functions, meaning when a node is about to be drawn, it is first sent to a layout function

# Set arguments to accept from the commandline
my_parser = argparse.ArgumentParser(description='This script is used to map a region.matrix file to a tyr phylogeny')

my_parser.add_argument('-f',
                       '--family',
                       metavar='FILE',
                       type=str,
                       help='a 2 column tsv with: tyrID, famID (optional)')

my_parser.add_argument('-t',
                       '--tree',
                       metavar='FILE',
                       type=str,
                       help='a newick tree')

my_parser.add_argument('-p',
                       '--prune',
                       metavar='FILE',
                       type=str,
                       help='a single column file of IDs to use to prune tree (optional)')

my_parser.add_argument('-a',
                       '--alignment',
                       metavar='FILE',
                       type=str,
                       help='a .fasta alignment file for mapping to the tree')

my_parser.add_argument('-s',
                       '--support',
                       metavar='STR',
                       type=str,
                       help='branch support in tree, either "iqtree" or "bootstrap"')


# Retrieve arguments from commandline
args = my_parser.parse_args()
family_file = args.family
tree_file = args.tree
alignment_file = args.alignment
prune_file = args.prune
supportType = args.support

# Check arguments exist
if tree_file is not None:
	if not os.path.isfile(tree_file):
		print('\nerror: the specified --tree file does not exist\n')
		sys.exit()
else:
	print('\nerror: please provide a --tree file\n')
	sys.exit()

if alignment_file is not None:
	if not os.path.isfile(alignment_file):
		print('\nerror: the specified --alignment_file file does not exist\n')
		sys.exit()
else:
	print('\nerror: please provide a --alignment_file file\n')
	sys.exit()

if supportType is not None:
	if not supportType == 'iqtree' and not supportType == 'bootstrap':
		print('\nerror: please specify either "iqtree" or "bootstrap" to --support \n')
		sys.exit()
else:
	print('\nerror: please provide a value to --support\n')
	sys.exit()


# Parse input family file
tyr2fam = {}
if family_file is not None:
	if not os.path.isfile(family_file):
		print('\nerror: the specified --family file does not exist\n')
		sys.exit()
	else:
		with open(family_file, 'r') as f:
			for line in f:
				splitLine = re.split(r'\t', line)
				tyr2fam[str(splitLine[0])] = str(splitLine[1].strip('\n'))


# Store input tree in TreeNode class, according to its specific newick format
# TIP: use format 1 for IQtree dual SH-alrt and rapid bootstrap values
# TIP: use format 2 for traditional bootstrap values
if supportType == 'bootstrap':
	tree = Tree(tree_file, format= 2) 
elif supportType == 'iqtree':
	tree = Tree(tree_file, format= 1) 

# Root tree at midpoint
R = tree.get_midpoint_outgroup()
tree.set_outgroup(R)

# Do an initial traverse of tree to collapse poorly supported branches
# TIP: from IQtree manual: one would typically start to rely on the clade if its SH-aLRT >= 80% and UFboot >= 95%
# TIP: traditional practice is to collapse all branches < 70% bootstrap support
for node in tree.traverse('preorder'):
	if not node.is_leaf():
		# For traditional bootstrap:
		if supportType == 'bootstrap':
			if node.support < 70:
				node.delete()
		elif supportType == 'iqtree':
			# For IQtree SH-aLRT and UFboot:
			if re.search('\/', node.name) and not node.is_leaf():
				splitsupport = re.split(r'\/', node.name)
				if float(splitsupport[0]) < 80 and float(splitsupport[1]) < 95:
					node.delete()

# prune tree, if a pruning file was provided
if args.prune is not None:

	# read in pruning file
	pruningIDs = []
	with open(prune_file, 'r') as f:
		for line in f:
			pruningIDs.append(line.strip('\n'))
		
	# prune tree
	tree.prune(pruningIDs, preserve_branch_length = True)

# creat a phylo-tree
treeTrimmed = tree.write(format=5)
treeAlignment = PhyloTree(treeTrimmed)
treeAlignment.link_to_alignment(alignment=alignment_file, alg_format="fasta")

# Create custom layout function, setting custom node styles and node faces
def mylayout(node):
	# If node is a leaf node
	if node.is_leaf():

		# Sets the base style of leaf nodes
		node.img_style["size"] = 0
		node.img_style["shape"] = "square"
		node.img_style["vt_line_width"] = 8
		node.img_style["hz_line_width"] = 8
		node.img_style["fgcolor"] = "black"
		leafcolor = 'black'

		# add MSA
		seqFace = faces.SeqMotifFace(seq=node.sequence, seqtype = 'aa', seq_format = 'seq', scale_factor=50, height=24, width=24)
		faces.add_face_to_node(seqFace, node, column=0, aligned=True)

		# add descriptions, if available
		if node.name in tyr2fam.keys():
			descFace = faces.TextFace(tyr2fam[node.name], fsize=24)
			descFace.margin_top = 1
			descFace.margin_bottom = 1
			descFace.border.margin = 1

			# Note that this faces is added in "aligned" mode
			faces.add_face_to_node(descFace, node, column=1, aligned=True)

		# Modify tip labels
		BigNameFace = faces.TextFace(node.name, fsize=24, fgcolor=leafcolor, ftype = 'Helvetica')
		faces.add_face_to_node(BigNameFace, node, column=0)
		

	# If node is an internal node
	else:
		# Sets the base style of internal nodes
		node.img_style["size"] = 0 
		node.img_style["shape"] = "square"
		node.img_style["vt_line_width"] = 8
		node.img_style["hz_line_width"] = 8
		node.img_style["fgcolor"] = "black"

		# add branch support values
		# for traditional bootstrap
		#support=faces.TextFace(int("{0:.0f}".format(node.support)),fsize=14, fgcolor="red", ftype = 'Helvetica')
		#faces.add_face_to_node(support, node, column=0, position="branch-top")

		# for IQTREE sh-ALRT and rapid bootstrap
		#support=faces.TextFace(node.name,fsize=6, fgcolor="black", ftype = 'Helvetica')
		#faces.add_face_to_node(support, node, column=0,position="branch-top")

# Customize tree style
treeAlignment.ladderize()

ts = TreeStyle()
ts.scale = 5000
ts.show_branch_support = False
ts.show_branch_length = False
ts.show_leaf_name = False
ts.draw_guiding_lines=True
ts.guiding_lines_type = 2
ts.guiding_lines_color = 'lightgray'

# Link layout function to tree style 
ts.layout_fn = mylayout


# Finally, open interactive viewer rendering with our tree style
treeAlignment.show(tree_style=ts)
#treeAlignment.render(tree_file+'.pdf', w = 12, h = 24, units = 'in', tree_style=ts)
