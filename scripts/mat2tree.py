#!/usr/bin/env python
import argparse
import os
import sys
import re
import colorsys
from ete3 import Tree, ClusterTree, faces, TreeStyle, RectFace, TextFace, add_face_to_node, AttrFace

# Adapted from http://etetoolkit.org/docs/latest/tutorial/tutorial_drawing.html#customizing-the-aspect-of-trees
# 	Image customization is performed through 4 main elements:
# 		1. Tree style:  controls the general aspect of the tree image
#		2. Node style:  controls various aspects of each single node
#		3. Node face:   graphical information that can be linked to nodes
#		4. Layout func: act as pre-drawing hooking functions, meaning when a node is about to be drawn, it is first sent to a layout function

# Set arguments to accept from the commandline
my_parser = argparse.ArgumentParser(description='This script is used to map a .mat file to a tyr phylogeny')

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

my_parser.add_argument('-m',
                       '--matrix',
                       metavar='FILE',
                       type=str,
                       help='a .mat file for mapping to the tree. First column header must be #Names')

my_parser.add_argument('-s',
                       '--support',
                       metavar='STR',
                       type=str,
                       help='branch support in tree, either "iqtree" or "bootstrap"')

my_parser.add_argument('-p',
                       '--prune',
                       metavar='FILE',
                       type=str,
                       help='a list of IDs in first column to use to prune tree (optional)')

my_parser.add_argument('-r',
                       '--root',
                       metavar='STR',
                       type=str,
                       help='a comma-delimited list of tips whose LCA will serve as root (optional)')


# Retrieve arguments from commandline
args = my_parser.parse_args()
family_file = args.family
tree_file = args.tree
matrix_file = args.matrix
supportType = args.support
prune_file = args.prune
rootTips = args.root

# Check arguments exist
if tree_file is not None:
	if not os.path.isfile(tree_file):
		print('\nerror: the specified --tree file does not exist\n')
		sys.exit()
else:
	print('\nerror: please provide a --tree file\n')
	sys.exit()

if matrix_file is not None:
	if not os.path.isfile(matrix_file):
		print('\nerror: the specified --matrix file does not exist\n')
		sys.exit()
else:
	print('\nerror: please provide a --matrix file\n')
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

header = []
with open(matrix_file, 'r') as f:
    first_line = f.readline()
    header = re.split(r'\t', first_line)
    header.pop(0)

# Store input tree in TreeNode class, according to its specific newick format
# TIP: use format 1 for IQtree dual SH-alrt and rapid bootstrap values
# TIP: use format 2 for traditional bootstrap values
if supportType == 'bootstrap':
	tree = Tree(tree_file, format= 2) 
elif supportType == 'iqtree':
	tree = Tree(tree_file, format= 1) 

# Root tree at midpoint
R = tree.get_midpoint_outgroup()
if rootTips is not None:
	focalTips = re.split(r',', rootTips)
	R = tree.get_common_ancestor(focalTips)
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
			splitLine = re.split(r'\t', line)
			pruningIDs.append(splitLine[0].strip('\n'))
		
	# prune tree
	tree.prune(pruningIDs, preserve_branch_length = True)

# remove support values (ClusterTree does not accept them)
# creat a matrix-tree
treeTrimmed = tree.write(format=5)
treeMatrix = ClusterTree(treeTrimmed, text_array=matrix_file)

# Create custom layout function, setting custom node styles and node faces
def mylayout(node):
	# If node is a leaf node
	if node.is_leaf():

		# Sets the base style of leaf nodes
		node.img_style["size"] = 0
		node.img_style["shape"] = "square"
		node.img_style["vt_line_width"] = 12
		node.img_style["hz_line_width"] = 12
		node.img_style["fgcolor"] = "black"
		leafcolor = 'black'

		# add descriptions, if available
		if node.name in tyr2fam.keys():
			descFace = faces.TextFace(tyr2fam[node.name], fsize=36)
			descFace.margin_top = 10
			descFace.margin_bottom = 10
			descFace.border.margin = 1

			# Note that this faces is added in "aligned" mode
			faces.add_face_to_node(descFace, node, column=0, aligned=True)

		# Modify tip labels
		BigNameFace = faces.TextFace(node.name, fsize=36, fgcolor=leafcolor, ftype = 'Helvetica')
		faces.add_face_to_node(BigNameFace, node, column=0)

	# If node is an internal node
	else:
		# Sets the base style of internal nodes
		node.img_style["size"] = 0 
		node.img_style["shape"] = "square"
		node.img_style["vt_line_width"] = 12
		node.img_style["hz_line_width"] = 12
		node.img_style["fgcolor"] = "black"

		# add branch support values
		# for traditional bootstrap
		#support=faces.TextFace(int("{0:.0f}".format(node.support)),fsize=14, fgcolor="red", ftype = 'Helvetica')
		#faces.add_face_to_node(support, node, column=0, position="branch-top")

		# for IQTREE sh-ALRT and rapid bootstrap
		#support=faces.TextFace(node.name,fsize=6, fgcolor="black", ftype = 'Helvetica')
		#faces.add_face_to_node(support, node, column=0,position="branch-top")



def setup_heatmap(tree, tree_style, header, center_value, color_up, color_down, color_center):
	DEFAULT_COLOR_SATURATION = 0.3
	BASE_LIGHTNESS = 1.1 #CHANGE here for brightness of heatmap, 2 is good for count matrices with larger max - min range
	def gradient_color(value, max_value, saturation=0.3, hue=0.1):    
		def rgb2hex(rgb):
			return '#%02x%02x%02x' % rgb
		def hls2hex(h, l, s):
			return rgb2hex( tuple(map(lambda x: int(x*255), colorsys.hls_to_rgb(h, l, s))))

		lightness = 1 - (value * BASE_LIGHTNESS) / max_value
		return hls2hex(hue, lightness, DEFAULT_COLOR_SATURATION)

    # Calculate max gradient value from the ClusterTree matrix
	maxv = abs(center_value - treeMatrix.arraytable._matrix_max)
	minv = abs(center_value - treeMatrix.arraytable._matrix_min)
	if center_value <= treeMatrix.arraytable._matrix_min:
		MAX_VALUE = minv + maxv
	else:
		MAX_VALUE = max(maxv, minv)
        
    # Add heatmap colors to tree
	for lf in tree:
		for i, value in enumerate(getattr(lf, "profile", [])):
			if value > center_value:
				color = gradient_color(abs(center_value - value), MAX_VALUE, hue=color_up)
			elif value < center_value:
				color = gradient_color(abs(center_value - value), MAX_VALUE, hue=color_down)
			else:
				color = color_center
			lf.add_face(RectFace(60, 60, color, color), position="aligned", column=i+1) #change the size of the heatmap here (width, height); column adds 1 if mapping a description to leaf
			# Uncomment to add numeric values to the matrix
			#lf.add_face(TextFace("%0.2f"%value, fsize=5), position="aligned", column=i)
# Add header 
	for i, name in enumerate(header):
		nameF = TextFace(name, ftype = "Helvetica", fsize=50)
		nameF.rotation = -90
		tree_style.aligned_header.add_face(nameF, column=i+1)

# This function adds faces and header labels to any ClusterTree, supporting
# heatmaps of 1 or 2 color series.  If center value is set to the minimum value
# observed in the matrix, only one color will be used used, otherwise values
# above the center will use a custom color gradient, and values bellow a
# different color gradient.  The base color of each gradient is specified with
# color_up and color_down variables, whose value should be between 0 and 1 (hue
# color value, see http://en.wikipedia.org/wiki/HSL_and_HSV)
		
# Customize tree style
treeMatrix.ladderize()

ts = TreeStyle()
ts.scale = 2000
ts.show_branch_support = False
ts.show_branch_length = False
ts.show_leaf_name = False
ts.draw_guiding_lines=True
ts.guiding_lines_type = 2
ts.guiding_lines_color = 'lightgray'

# Link layout function to tree style 
ts.layout_fn = mylayout

# heatmap setup
setup_heatmap(treeMatrix, ts, header, center_value=0, color_up=0.5, color_down=0.2, color_center="ghostwhite")

# Finally, open interactive viewer rendering with our tree style
treeMatrix.render(tree_file+'.pdf', w = 12, h = 26, units = 'in', tree_style=ts)
treeMatrix.show(tree_style=ts)
