#!/usr/bin/env python3

from ete3 import Tree, TreeStyle, NodeStyle, faces, PieChartFace, TextFace, CircleFace
import sys
import math
from operator import add
from matplotlib import colors as mpl_colors
from matplotlib import cm as mpl_colormap
from datetime import datetime


class TreeMaker(object):

	def __init__(self, tree, treefile, prune, sqrt = False, collapse = True, show_empty_species = True, startnode = None, countTotal = False):

		if startnode is None:
			self.t = tree
		else:
			node = tree.search_nodes(name=startnode)
			if not node:
				raise KeyError('Node {} not found in tree.'.format(startnode))
			self.t = node[0]

		self.sqrt = sqrt
		self.collapse = collapse
		self.empty = show_empty_species
		self.treefile = treefile
		self.prune = prune
		self.startnode = startnode
		self.count = countTotal

		## Colors for the features
		##                         0          1          2          3          4          5          6          7
		##                         A          B          C         A+B        A+C        B+C        all        none
		self.standardcolors = ['#AA0000', '#00AA00', '#55AAFF', '#DDDD00', '#AA00AA', '#3333FF', '#000000', '#c8c8c8']

		## some more colors from colorbrewer2.org, last set edited a bit
		## This allows a better of different colors & additional support for the empty feature
		self.emptycolor = '#c8c8c8'
		self.colors_small = ['#e41a1c','#377eb8','#4daf4a','#984ea3','#ff7f00','#ffff33']
		self.colors_big = ['#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99','#e31a1c','#fdbf6f','#ff7f00','#cab2d6','#6a3d9a','#ffffb3', '#ffea4d', '#b15928', '#000000']

		## Data should be read into self.elements (dict)

		## These functions should be called (in order) at the end of the init() of each *child* class
		## Optional TODO: always call this init *after* child specific init part and thereby run the following functions

		#self.addFeatures()
		#self.readPruningFile()
		#self.pruneTree()
		#self.treeLayout()

		## This code is helpful for if using an arbitrary number of features:

		# # use different range of colors depeding on number of features
		# if len(self.featurelist) <= 6:
		# 	self.colors = self.colors_small
		# elif len(self.featurelist) <= 14:
		# 	self.colors = self.colors_big
		# else:
		#	self.colors = self.colors_big
		# 	print('Number of features exceeds maximum number of supported colors (14)! Reusing colors.')
		# 	while len(self.colors) < len(self.featurelist):
		# 		self.colors += self.colors_big
		# # make self.colors as long as the featurelist + grey to support empty
		# self.colors = self.colors[0:len(self.featurelist)]
		# self.colors += [self.emptycolor]

	# end __init__()


	# Exemplary minimal function to be subclassed, (better) keep the non-commented code from here
	# Optional TODO: could be generalised by using the self.featurelist, thereby restricting all subclasses to PieCharts
	def layout(self, node):
		######  How to make PieCharts:	############################################
		# try:
		# 	percents = [round(100.0*node.dummyfeature/node.total, round(100.0*node.f_none/node.total)]
		# except ZeroDivisionError:
		# 	txt = TextFace(' ' + node.plainName, ftype='Arial', fgcolor='#000000')
		# else:
		# 	i=0
		# 	if sum(percents) != 100:
		# 		x = sorted(percents)[-1]
		# 		i = percents.index(x)
		# 		percents[i] += 100 - sum(percents)
        #
		# 	if self.sqrt:
		# 		size = int(math.sqrt(node.total)*5+1)
		# 	else:
		# 		size = int(math.log10(node.total)*20+10)
        #
		# 	P = PieChartFace(percents, size, size, self.colors)
		# 	faces.add_face_to_node(P, node, 0, position="branch-top")
        #
		fgcol = '#000000'

		if self.count:
			cnt = node.total
			txt = TextFace(' {}\n {:d}'.format(node.plainName, cnt), ftype='Arial', fgcolor=fgcol)
		else:
			txt = TextFace(' ' + node.plainName, ftype='Arial', fgcolor=fgcol)
		#############################################################################

		if node.is_leaf():
			faces.add_face_to_node(txt, node, 0, position="branch-right")
		else:
			faces.add_face_to_node(txt, node, 0, position="branch-bottom")

	# end layout


	# Exemplary minimal function to be subclassed, (better) keep the non-commented code from here
	# Optional TODO: could be generalised by using the self.featurelist for the legend and adding a list for the used files
	def treeLayout(self):
		# Add a tree style and a legend

		self.ts = TreeStyle()
		self.ts.show_scale = False
		self.ts.show_leaf_name = False
		self.ts.extra_branch_line_type = 0
		self.ts.extra_branch_line_color = 'black'
		self.ts.layout_fn = self.layout
		self.ts.branch_vertical_margin = -10

		self.ts.legend_position = 4

		self.ts.legend.add_face(TextFace(datetime.now().strftime('%a, %d.%m.%Y; %H:%M:%S'), fsize=8), column=1)
		self.ts.legend.add_face(TextFace('Tree: ' + self.treefile, fsize=8), column=1)
		self.ts.legend.add_face(TextFace('Pruning: ' + self.prune, fsize=8), column=1)

		# Order the graphical tree output from highly branched nodes on top and less branched nodes to the bottom
		self.t.ladderize(1)

	# end treeLayout


	# Exemplary minimal function to be subclassed, (better) keep the non-commented code from here
	def addFeatures(self):
		# Traverse the tree and add features to each node
		for n in self.t.traverse('postorder'):
			if n.name == 'NoName' or not n.name:
				n.name = 'Biota^1'

			# This requires that the Newick files nodes named: 'OraganimName^Taxid'
			# TODO: add excpetion if Newick file has different format
			n.add_features(plainName = n.name.split('^')[0])
			n.add_features(taxid = int(n.name.split('^')[1]))

			#Initialise desired features
			#n.add_features(**kwargs)
			n.add_features(dummy = 0)
			n.add_features(f_none = 0)
			n.add_features(total = 0)

			if n.is_leaf():
				countThis = True

				# Here the actual data has to be loaded (best: self.elements) and used to set features depending on TaxID
				# The features names used on the nodes should/can also be present on the base class (self).
				# This way the loaded data or the self.features can be containers for the taxIDs or data-values to be used/checked

				#Replace this part
				dummydata = set()
				if dummydata:
					n.dummy += 1
				#Keep this part
				elif self.empty:
					n.f_none += 1
				else:
					countThis = False

				if countThis:
					n.total += 1

			else:
				for x in n.children:
					#Add all your features here for propagation
					n.dummy += x.dummy
					#Always keep these two
					n.f_none += x.f_none
					n.total += x.total

	# end addFeatures()


	# Universal
	def readPruningFile(self):
		self.pruneAfter = set()
		self.pruneAfterChildren = set()
		self.dontPrune = set()
		self.detach = set()
		self.delete = set()

		try:
			with open(self.prune, 'r') as f:
				active = self.pruneAfter
				for line in f:
					line = line.strip()

					if not line:
						continue

					if line.startswith('#'):
						continue

					if line.startswith('!'):
						keyword = line[1:].lower()
						if keyword == 'pruneafter':
							active = self.pruneAfter
						elif keyword == 'pruneafterchildren':
							active = self.pruneAfterChildren
						elif keyword == 'dontprune':
							active = self.dontPrune
						elif keyword == 'detach':
							active = self.detach
						elif keyword == 'delete':
							active = self.delete
						else:
							print('Keyword {} not found!'.format(keyword))

						continue

					line = line.split('#')[0].rstrip()

					if '^' in line:
						active.add(int(line.split('^')[1]))
					else:
						try:
							active.add(int(line))
						except ValueError:
							active.add(line)

		except IOError:
			print('File {} not found!'.format(self.prune))

	# end readPruningFile()


	# Universal
	def pruneTree(self):
		# Start pruning the tree

		# First, we need to collect the nodes to prune and afterwards prune them.
		# Otherwise, tree traversing doesn't work.
		toPrune = []
		toDetach = []
		toFinalDelete = []
		for n in self.t.traverse('postorder'):
			# Prune after children
			if n.taxid in self.pruneAfterChildren or n.plainName in self.pruneAfterChildren:
				for m in n.children:
					toPrune.append(m)

			# Prune directly
			if n.taxid in self.pruneAfter or n.plainName in self.pruneAfter:
				toPrune.append(n)

			# Detach directly
			if n.taxid in self.detach or n.plainName in self.detach:
				toDetach.append(n)

			# Delete in the end
			if n.taxid in self.delete or n.plainName in self.delete:
				toFinalDelete.append(n)

		# Now, we have the nodes to keep. We need a list of nodes to detach.
		for n in toPrune:
			for m in n.children:
				toDetach.append(m)

		# Finally, we can detach the nodes.
		for n in toDetach:
			if n.taxid not in self.dontPrune and n.plainName not in self.dontPrune:
				n.detach()

		# Now that the tree is pruned, we may go through it ones again and see, whether
		# we can collapse some nodes (if wanted by user)
		if self.collapse:
			toDelete = []
			for n in self.t.traverse('postorder'):

				# Node attributes that should lead to collapsing (= no information in node) [marked *]
				#
				# * len(n.children) == 1 -> reduce long singular lineages to one step, leaves dont have children
				# 	- if not self.empty only children with n.total > 0 have to be considered: len([x for x in n.children if x.total])
				#	- the root node can be 'deleted', therefore at that point complete collapsing fails right now
				#
				# * n.total == 0 [and n not in toPrune] -> empty nodes, (Pruning makes new leaves, keeping these) [^1]
				#	- this obviously doesn't apply showing empty nodes
				#	- I think in principle pruning commands should not affect collapsing (collapse should only make tree more
				#     readable, but dont change information). If (default) empty nodes should disappear, what about nodes that
				#     are empty & pruned ? Keep them just because they were mentioned ? (currently for leaves made by pruning ^1).
				#     It would probably be more consequent to collapse emtpy nodes regardless of pruning.
				#	  Maybe a '!dontcollapse' option for specific nodes could be introduced for the pruningfile


				# not sure what exactly the originally included 'len(n.get_sisters())==0' was supposed to do, but it caused too much collapsing

				if n.is_leaf():
					if not self.empty and n.total == 0 and n not in toPrune:
						toDelete.append(n)
				elif self.empty: # with and len(..) not sure if this could jump to the next elif if the part after and is False
					if len(n.children) == 1:
						toDelete.append(n)
				#not self.empty
				elif len([x for x in n.children if x.total]) == 1 or n.total == 0:
					toDelete.append(n)

			for n in toDelete:
				n.delete()

		for n in toFinalDelete:
			n.delete()

	# end pruneTree()


	# Universal
	def showTree(self):
		self.t.show(tree_style=self.ts)


	# Universal
	def renderTree(self, outfile, width, dpi):
		self.t.render(outfile, w = width, dpi = dpi, tree_style = self.ts)

# end TreeMaker base class


class Combinations(TreeMaker):
	def __init__(self, config, attr, **kwargs):
		TreeMaker.__init__(self, **kwargs)

		self.config = config
		self.attr = attr

		self.colors = self.standardcolors

		# Read tree attributes into dict
		self.elements = {}
		with open(self.attr, 'r') as f:
			for line in f:
				l = line.rstrip().split('\t')
				self.elements[int(l[0])] = set([l[1][i:i + 2] for i in range(0, len(l[1]), 2)])
			#      ^ `elements` looks like this: {1234: {'aA', 'aF', 'bE'}, ...}

		self.f_a = None  # Code for feature 1
		self.f_b = None  # Code for feature 2
		self.f_c = None  # Code for feature 3

		self.name_a = None  # Name for feature 1
		self.name_b = None  # Name for feature 2
		self.name_c = None  # Name for feature 3

		self.counterElements = []


		self.readConfig()

		self.addFeatures()

		self.readPruningFile()

		self.pruneTree()

		self.treeLayout()

# for config/feature-combinations
	def readConfig(self):
		with open(self.config, 'r') as f:
			for line in f:
				line = line.strip()

				if line.startswith('!'):
					line = line[1:].replace(',', ' ')
					self.colors = line.split()
					l = len(self.colors)
					if l > 8:
						print('Too many colors in config file')
						self.colors = self.colors[:8]
					elif l < 8:
						print('Too few colors in config file')
						self.colors.extend(self.standardcolors[l:])
					continue

				line = line.split('#')[0].strip()

				if not line:
					continue

				if line.startswith('>'):
					line = line[1:]
					self.counterElements = [line[i:i+2] for i in range(0, len(line), 2)]
					continue

				line = line.split(':')
				if self.f_a is None:
					self.f_a = line[0]
					self.name_a = line[1]
					continue

				if self.f_b is None:
					self.f_b = line[0]
					self.name_b = line[1]
					continue

				if self.f_c is None:
					self.f_c = line[0]
					self.name_c = line[1]
					continue

	# end readConfig()


	def layout(self, node):
		try:
			percents = [round(100.0*node.f_a/node.total), round(100.0*node.f_b/node.total), round(100.0*node.f_c/node.total), round(100.0*node.f_ab/node.total), round(100.0*node.f_ac/node.total), round(100.0*node.f_bc/node.total), round(100.0*node.f_all/node.total), round(100.0*node.f_none/node.total)]
		except ZeroDivisionError:
			txt = TextFace(' ' + node.plainName, ftype='Arial', fgcolor='#000000')
		else:
			i=0
			if sum(percents) != 100:
				x = sorted(percents)[-1]
				i = percents.index(x)
				percents[i] += 100 - sum(percents)

			if self.sqrt:
				size = int(math.sqrt(node.total)*5+1)
			else:
				size = int(math.log10(node.total)*20+10)

			P = PieChartFace(percents, size, size, self.colors)
			faces.add_face_to_node(P, node, 0, position="branch-top")

			#fgcol = mpl_colors.rgb2hex(mpl_colormap.brg(col/10.0)) # foreground/text color depending on NAD consumer density
			fgcol = '#000000'

			if node.counter:
				cnt = float(sum(node.counter))/float(node.total)
				txt = TextFace(' {}\n {:.2f}'.format(node.plainName, cnt), ftype='Arial', fgcolor=fgcol)
			elif self.count:
				cnt = node.total
				txt = TextFace(' {}\n {:d}'.format(node.plainName, cnt), ftype='Arial', fgcolor=fgcol)
			else:
				txt = TextFace(' ' + node.plainName, ftype='Arial', fgcolor=fgcol)
			#txt.margin_right = -30

		if node.is_leaf():
			faces.add_face_to_node(txt, node, 0, position="branch-right")
		else:
			faces.add_face_to_node(txt, node, 0, position="branch-bottom")

	# end layout


	def treeLayout(self):
		# Add a tree style and a legend

		self.ts = TreeStyle()
		self.ts.show_scale = False
		self.ts.show_leaf_name = False
		self.ts.extra_branch_line_type = 0
		self.ts.extra_branch_line_color = 'black'
		self.ts.layout_fn = self.layout
		self.ts.branch_vertical_margin = -10

		self.ts.legend_position = 4

		if self.name_a:
			self.ts.legend.add_face(CircleFace(10, self.colors[0]), column=0)
			self.ts.legend.add_face(TextFace(self.name_a), column=1)

		if self.name_b:
			self.ts.legend.add_face(CircleFace(10, self.colors[1]), column=0)
			self.ts.legend.add_face(TextFace(self.name_b), column=1)

		if self.name_c:
			self.ts.legend.add_face(CircleFace(10, self.colors[2]), column=0)
			self.ts.legend.add_face(TextFace(self.name_c), column=1)

		if self.name_a and self.name_b:
			self.ts.legend.add_face(CircleFace(10, self.colors[3]), column=0)
			self.ts.legend.add_face(TextFace(self.name_a + ' + ' + self.name_b), column=1)

		if self.name_a and self.name_c:
			self.ts.legend.add_face(CircleFace(10, self.colors[4]), column=0)
			self.ts.legend.add_face(TextFace(self.name_a + ' + ' + self.name_c), column=1)

		if self.name_b and self.name_c:
			self.ts.legend.add_face(CircleFace(10, self.colors[5]), column=0)
			self.ts.legend.add_face(TextFace(self.name_b + ' + ' + self.name_c), column=1)

		if self.name_a and self.name_b and self.name_c:
			self.ts.legend.add_face(CircleFace(10, self.colors[6]), column=0)
			self.ts.legend.add_face(TextFace(self.name_a + ' + ' + self.name_b + ' + ' + self.name_c), column=1)

		if self.empty:
			self.ts.legend.add_face(CircleFace(10, self.colors[7]), column=0)
			self.ts.legend.add_face(TextFace('none'), column=1)

		self.ts.legend.add_face(TextFace(datetime.now().strftime('%a, %d.%m.%Y; %H:%M:%S'), fsize=8), column=1)
		self.ts.legend.add_face(TextFace('Tree: ' + self.treefile, fsize=8), column=1)
		self.ts.legend.add_face(TextFace('Attributes: ' + self.attr, fsize=8), column=1)
		self.ts.legend.add_face(TextFace('Config: ' + self.config, fsize=8), column=1)
		self.ts.legend.add_face(TextFace('Pruning: ' + self.prune, fsize=8), column=1)

		# Order the graphical tree output from highly branched nodes on top and less branched nodes to the bottom
		self.t.ladderize(1)

	# end treeLayout


	def addFeatures(self):
		# Traverse the tree and add features to each node
		for n in self.t.traverse('postorder'):
			if n.name == 'NoName' or not n.name:
				n.name = 'Biota^1'

			n.add_features(plainName = n.name.split('^')[0])
			n.add_features(taxid = int(n.name.split('^')[1]))

			n.add_features(f_a = 0)
			n.add_features(f_b = 0)
			n.add_features(f_c = 0)
			n.add_features(f_ab = 0)
			n.add_features(f_ac = 0)
			n.add_features(f_bc = 0)
			n.add_features(f_all = 0)
			n.add_features(f_none = 0)
			n.add_features(total = 0)
			n.add_features(counter = [0]*len(self.counterElements))
			if n.is_leaf():
				countThis = True
				s = self.elements[n.taxid]

				if self.f_a in s and self.f_b in s and self.f_c in s:
					n.f_all += 1
				elif self.f_b in s and self.f_c in s:
					n.f_bc += 1
				elif self.f_a in s and self.f_c in s:
					n.f_ac += 1
				elif self.f_a in s and self.f_b in s:
					n.f_ab += 1
				elif self.f_c in s:
					n.f_c += 1
				elif self.f_b in s:
					n.f_b += 1
				elif self.f_a in s:
					n.f_a += 1
				elif self.empty:
					n.f_none += 1
				else:
					countThis = False

				if countThis:
					n.total += 1

					for i, c in enumerate(self.counterElements):
						if c in s:
							n.counter[i] += 1

			else:
				for x in n.children:
					n.f_a += x.f_a
					n.f_b += x.f_b
					n.f_c += x.f_c
					n.f_ab += x.f_ab
					n.f_ac += x.f_ac
					n.f_bc += x.f_bc
					n.f_all += x.f_all
					n.f_none += x.f_none
					n.total += x.total

					n.counter = list(map(add, n.counter, x.counter))

	# end addFeatures()

# end Combinations class


## While originally intended for crosshits, this basically takes protein IDs + taxIDs with e-values & shows the distribution. Rename it ??
class Crosshits(TreeMaker):
	def __init__(self, combo, bins, **kwargs):
		TreeMaker.__init__(self, **kwargs)

		self.combo = combo

		# Read tree attributes into dict
		self.elements = {}
		with open(self.combo, 'r') as f:
			for line in f:
				line = line.rstrip('\n').split('#')[0]
				if not line:
					continue

				ev, num, entries = line.split('\t')

				for entry in entries.split(','):
					if not entry:
						continue
					acc = entry.split('^')[0]
					tax = int(entry.split('^')[1])
					if tax in self.elements:
						self.elements[tax] += [int(ev)]
					else:
						self.elements[tax] = [int(ev)]

		# self.featurelist list strs with feature names to be used (f_0, f_1, ...)
		# the features of the top-level class (self.f_0, ...) equal a set-range of corresponding evalues (exponents only)

		# add Support for complete custum e-value ranges:
		if isinstance(bins, set):
			self.featurelist = ['f_{:02d}'.format(i) for i in range(len(bins) - 1)]
			borders = sorted(list(bins))
			for i, f in enumerate(self.featurelist):
				setattr(self, f, set(range(borders[i], borders[i + 1])))
			# add last border to last feature-set
			lastf = getattr(self, self.featurelist[-1])
			lastf.add(borders[-1])

		# autocalcualate bins for evalue ranges
		else:
			self.featurelist = ['f_{:02d}'.format(i) for i in range(bins[2])]

			steps = int((bins[1] - bins[0]) / bins[2])
			for i, f in enumerate(self.featurelist):
				setattr(self, f, set(range(bins[0] + steps * i, bins[0] + steps * (i + 1))))
			# extend last range up to bins[1], more than last one may be missing due to rounding
			lastf = getattr(self, self.featurelist[-1])
			lastf.update(range(bins[1] - steps, bins[1] + 1))

		#use different range of colors depeding on number of features
		if len(self.featurelist) <= 6:
			self.colors = self.colors_small
		elif len(self.featurelist) <= 14:
			self.colors = self.colors_big
		else:
			print('Number of features exceeds maximum number of supported colors (14)! Reusing colors.')
			self.colors = self.colors_big
			while len(self.colors) < len(self.featurelist):
				self.colors += self.colors_big
		#make self.colors as long as the featurelist + grey to support empty
		self.colors = self.colors[0:len(self.featurelist)]
		self.colors += [self.emptycolor]

		self.featurenames = ['e-values between 1e-{} and 1e-{}'.format(min(getattr(self, f)), max(getattr(self, f))) for f in self.featurelist]

		self.addFeatures()

		self.readPruningFile()

		self.pruneTree()

		self.treeLayout()

	# end __init__()


	def layout(self, node):
		try:
			percents = [round(100.0 * getattr(node, f) / node.total) for f in self.featurelist] + [round(100.0 * node.f_none / node.total)]
		except ZeroDivisionError:
			txt = TextFace(' ' + node.plainName, ftype='Arial', fgcolor='#000000')
		else:
			if sum(percents) != 100:
				x = sorted(percents)[-1]
				i = percents.index(x)
				percents[i] += 100 - sum(percents)

			if self.sqrt:
				size = int(math.sqrt(node.total) * 5 + 1)
			else:
				size = int(math.log10(node.total) * 20 + 10)

			P = PieChartFace(percents, size, size, self.colors)
			faces.add_face_to_node(P, node, 0, position="branch-top")

			fgcol = '#000000'

			if self.count:
				cnt = node.total
				txt = TextFace(' {}\n {:d}'.format(node.plainName, cnt), ftype='Arial', fgcolor=fgcol)
			else:
				txt = TextFace(' ' + node.plainName, ftype='Arial', fgcolor=fgcol)

		if node.is_leaf():
			faces.add_face_to_node(txt, node, 0, position="branch-right")
		else:
			faces.add_face_to_node(txt, node, 0, position="branch-bottom")

	# end layout


	def treeLayout(self):
		# Add a tree style and a legend

		self.ts = TreeStyle()
		self.ts.show_scale = False
		self.ts.show_leaf_name = False
		self.ts.extra_branch_line_type = 0
		self.ts.extra_branch_line_color = 'black'
		self.ts.layout_fn = self.layout
		self.ts.branch_vertical_margin = -10

		self.ts.legend_position = 4

		for name, color in zip(self.featurenames, self.colors):
			self.ts.legend.add_face(CircleFace(10, color), column=0)
			self.ts.legend.add_face(TextFace(name), column=1)

		if self.empty:
			self.ts.legend.add_face(CircleFace(10, self.colors[-1]), column=0)
			self.ts.legend.add_face(TextFace('none'), column=1)

		self.ts.legend.add_face(TextFace(datetime.now().strftime('%a, %d.%m.%Y; %H:%M:%S'), fsize=8), column=1)
		self.ts.legend.add_face(TextFace('Tree: ' + self.treefile, fsize=8), column=1)
		self.ts.legend.add_face(TextFace('Crosshits: ' + self.combo[10:], fsize=8), column=1)
		self.ts.legend.add_face(TextFace('Pruning: ' + self.prune, fsize=8), column=1)

		# Order the graphical tree output from highly branched nodes on top and less branched nodes to the bottom
		self.t.ladderize(1)

	# end treeLayout


	def addFeatures(self):
		# Traverse the tree and add features to each node

		for n in self.t.traverse('postorder'):
			if n.name == 'NoName' or not n.name:
				n.name = 'Biota^1'

			n.add_features(plainName=n.name.split('^')[0])
			n.add_features(taxid=int(n.name.split('^')[1]))

			kwargs = {f: 0 for f in self.featurelist}
			n.add_features(**kwargs)

			n.add_features(f_none=0)
			n.add_features(total=0)



			if n.is_leaf():
				countThis = True
				i = 0
				try:
					# Tree is derived from one protein, therefore it has more entries than crosshits
					# There may be multiple evalues ('cross-crosshits') per species
					evs = self.elements[n.taxid].copy()

					while evs:
						ev = evs.pop()
						i += 1
						for f in self.featurelist:
							if ev in getattr(self, f):
								# n.'f' += 1
								setattr(n, f, getattr(n, f) + 1)
								break

				except KeyError:
					# no crosshit at this leaf
					if self.empty:
						i += 1
						n.f_none += 1
					else:
						countThis = False

				if countThis:
					n.total += i

			else:
				for x in n.children:
					for f in self.featurelist:
						# n.'f' += x.'f'
						setattr(n, f, getattr(n, f) + getattr(x, f))

					n.f_none += x.f_none
					n.total += x.total

	# end addFeatures()

# end Crosshits class


def greyscale(n):
	'''
	Function to make a list of greyscale colors
	:param n: number of colors wanted
	:return: list with colors in hex format
	'''
	vals = [255./n * i for i in range(n)]
	#convert to hex & add color for empty (bright red)
	colors = ['#{0:02x}{0:02x}{0:02x}'.format(int(i)) for i in vals[::-1]] + ['#fb9a99']
	return colors


# ToDo: add better description of file properties - where?
# supports comments(#), one cluster per line (tax or acc^tax, comma separated), !Clustername in extra line above the actual Cluster
class Clusters(TreeMaker):
	def __init__(self, filename, dropclusters, usegrey=False, **kwargs):
		TreeMaker.__init__(self, **kwargs)

		self.loadfile = filename
		self.dontuse = dropclusters

		self.colors = self.standardcolors

		# Read tree attributes (= clusters) from file
		# self.featurelist list strs with feature names to be used (f_00, f_01, ...)

		# the features of the top-level class (self.f_00, ...) poinst to a list with the accession ids for this cluster
		# this is actually not even needed - taxids in clusters are the only thing the script needs

		self.elements = {}
		self.featurelist = []
		self.featurenames = {}

		with open(self.loadfile, 'r') as f:
			i = 0
			for line in f:
				line = line.rstrip().split('#')[0]
				if not line:
					continue

				if i in self.dontuse:
					if not line.startswith('!'):
						i += 1
					continue

				#Allows definition of cluster names in file
				if line.startswith('!'):
					self.featurenames['f_{:02d}'.format(i)] = line[1:]
					continue

				cluster = line.split(',')
				for entry in cluster:
					if '^' in entry:
						tax = int(entry.split('^')[1])
					else:
						tax = int(entry)

					if tax in self.elements:
						self.elements[tax] += ['f_{:02d}'.format(i)]
					else:
						self.elements[tax] = ['f_{:02d}'.format(i)]

				self.featurelist.append('f_{:02d}'.format(i))
				setattr(self, 'f_{:02d}'.format(i), [n.split('^')[0] for n in cluster])

				i += 1

		if usegrey:
			self.colors = greyscale(len(self.featurelist))
		else:
			#use different range of colors depeding on number of features
			if len(self.featurelist) <= 6:
				self.colors = self.colors_small
			elif len(self.featurelist) <= 14:
				self.colors = self.colors_big
			else:
				print('Number of features exceeds maximum number of supported colors (14)! Reusing colors.')
				self.colors = self.colors_big
				while len(self.colors) < len(self.featurelist):
					self.colors += self.colors_big
			#make self.colors as long as the featurelist + add empty color
			self.colors = self.colors[0:len(self.featurelist)]
			self.colors += [self.emptycolor]

		for f in self.featurelist:
			if f not in self.featurenames:
				self.featurenames[f] ='Cluster Nr: {}'.format(f[2:])

		self.addFeatures()

		self.readPruningFile()

		self.pruneTree()

		self.treeLayout()

	# end __init__()


	def layout(self, node):
		try:
			percents = [round(100.0 * getattr(node, f) / node.total) for f in self.featurelist] + [round(100.0 * node.f_none / node.total)]
		except ZeroDivisionError:
			txt = TextFace(' ' + node.plainName, ftype='Arial', fgcolor='#000000')
		else:
			if sum(percents) != 100:
				x = sorted(percents)[-1]
				i = percents.index(x)
				percents[i] += 100 - sum(percents)

			if self.sqrt:
				size = int(math.sqrt(node.total) * 5 + 1)
			else:
				size = int(math.log10(node.total) * 20 + 10)

			P = PieChartFace(percents, size, size, self.colors)
			faces.add_face_to_node(P, node, 0, position="branch-top")

			fgcol = '#000000'

			if self.count:
				cnt = node.total
				txt = TextFace(' {}\n {:d}'.format(node.plainName, cnt), ftype='Arial', fgcolor=fgcol)
			else:
				txt = TextFace(' ' + node.plainName, ftype='Arial', fgcolor=fgcol)

		if node.is_leaf():
			faces.add_face_to_node(txt, node, 0, position="branch-right")
		else:
			faces.add_face_to_node(txt, node, 0, position="branch-bottom")

	# end layout


	def treeLayout(self):
		# Add a tree style and a legend

		self.ts = TreeStyle()
		self.ts.show_scale = False
		self.ts.show_leaf_name = False
		self.ts.extra_branch_line_type = 0
		self.ts.extra_branch_line_color = 'black'
		self.ts.layout_fn = self.layout
		self.ts.branch_vertical_margin = -10

		self.ts.legend_position = 4

		for feature, color in zip(sorted(self.featurenames.keys()), self.colors):
			self.ts.legend.add_face(CircleFace(10, color), column=0)
			self.ts.legend.add_face(TextFace(self.featurenames[feature]), column=1)

		if self.empty:
			self.ts.legend.add_face(CircleFace(10, self.colors[-1]), column=0)
			self.ts.legend.add_face(TextFace('none'), column=1)

		self.ts.legend.add_face(TextFace(datetime.now().strftime('%a, %d.%m.%Y; %H:%M:%S'), fsize=8), column=1)
		self.ts.legend.add_face(TextFace('Tree: ' + self.treefile, fsize=8), column=1)
		self.ts.legend.add_face(TextFace('Clusters: ' + self.loadfile, fsize=8), column=1)
		self.ts.legend.add_face(TextFace('Pruning: ' + self.prune, fsize=8), column=1)

		# Order the graphical tree output from highly branched nodes on top and less branched nodes to the bottom
		self.t.ladderize(1)

	# end treeLayout


	def addFeatures(self):
		# Traverse the tree and add features to each node
		for n in self.t.traverse('postorder'):
			if n.name == 'NoName' or not n.name:
				n.name = 'Biota^1'

			n.add_features(plainName=n.name.split('^')[0])
			n.add_features(taxid=int(n.name.split('^')[1]))

			kwargs = {f: 0 for f in self.featurelist}
			n.add_features(**kwargs)

			n.add_features(f_none=0)
			n.add_features(total=0)

			if n.is_leaf():
				countThis = True
				i = 0
				try:
					# Tree is derived from blast hits, number of promters is lower, self.elements doesnt have values for all tax-ids
					features = self.elements[n.taxid].copy()

					while features:
						f = features.pop()
						i += 1
						setattr(n, f, getattr(n, f) + 1)

				except KeyError:
					# No match in self.elements
					if self.empty:
						n.f_none += 1
						i += 1
					else:
						countThis = False

				if countThis:
					n.total += i

			else:
				for x in n.children:
					for f in self.featurelist:
						# n.'f' += x.'f'
						setattr(n, f, getattr(n, f) + getattr(x, f))

					n.f_none += x.f_none
					n.total += x.total

	# end addFeatures()

# end Clusters class


if __name__ == '__main__':
	import argparse

	parser = argparse.ArgumentParser(description='Create a tree with piecharts. Different modes are supported: *combination* of proteins, *crosshit* e-value distribution, Protein *cluster* distributons.')

	parser.add_argument('-m', '--mode', default='com', help='What to plot: combination, crosshits, or cluster. Default is combination.')

	# General options
	parser.add_argument('-t', '--tree', default='trees/general.tre', help='Filename of .tre file for the tree. The file must be in Newick format and the names in the format: name^taxid')
	parser.add_argument('-p', '--prune', default='tree_to_prune.txt', help='Filename of the file with information, which nodes to prune and delete.')
	parser.add_argument('--startnode', default=None, help='If the tree should not start from root, give the exact name of the node that should be the root of the tree.')
	parser.add_argument('--sqrt', action='store_true', help='Give this flag to use square root instead of log10 for the circle size')
	parser.add_argument('-n', '--nocollapse', action='store_true', help='If this flag is given, dont collapse the tree.')
	parser.add_argument('-e', '--empty', action='store_true', help='If this flag is given, show also species that do not have any of the given features.')
	parser.add_argument('-s', '--show', action='store_true', help='Show tree in ETE Tree Viewer (good for inspection)')
	parser.add_argument('-o', '--outfile', default='', help='If a filename is given, save the tree to that file. Valid extensions are: svg, pdf, png')
	parser.add_argument('-w', '--width', type=int, default=1000, help='Width of the resulting picture in pixels.')
	parser.add_argument('--count_total', action='store_true', help='Print the total number of hits in addition to Piechart size.')

	# Options only for 'combination'
	parser.add_argument('-c', '--config', default='tree_config.txt', help='Filename of the file that contains the config with the tree features and the colors. [combinations only]')
	parser.add_argument('-a', '--attr', default='attributes.txt', help='Filename of the file with the tree attributes. [combinations only]')

	# Options only for 'crosshit'
	parser.add_argument('--proteins', help='Filename of the .tsv file with the crosshits of two proteins (generated by phylogenetics.py step 12). [crosshits only]')
	parser.add_argument('-b', '--bins', type=int, default=[5], nargs=1, help='Number of bins/features to auto-creation of e-value ranges [crosshits only]')
	parser.add_argument('-r', '--ranges', type=int, nargs='+', default=[30,150], help='Min. 2 values. Either lower & upper bounds for ranges autocalculation, or multiple range borders [crosshits only]')

	# Options only for 'cluster'
	parser.add_argument('--clusters', help='Filename of the .csv file that contains one cluster per line (comma separated TaxID/acc^TaxID) [cluster only].')
	parser.add_argument('-d', '--drop', nargs='+', type=int, default=[], help='Indices of clusters that should not be included (First index = 0). [cluster only]')
	parser.add_argument('-g', '--greyscale', action='store_true', help='Use greyscale colors instead of normal ones. [cluster only]')


	args = parser.parse_args()

	t = Tree(args.tree, format=8)

	kwargs = {'tree': t, 'treefile': args.tree, 'prune': args.prune, 'sqrt': args.sqrt, 'collapse': not args.nocollapse, 'show_empty_species': args.empty, 'startnode': args.startnode, 'countTotal': args.count_total}

	if args.mode.startswith('com'):
		kwargs['config'] = args.config
		kwargs['attr'] = args.attr
		TM = Combinations(**kwargs)

	elif args.mode.startswith('cross'):
		if len(args.ranges) == 2:
			bins = args.ranges + args.bins
		elif len(args.ranges) > 2:
			bins = set(args.ranges)
		else:
			print('ranges option (-r) needs at least two values!', file=sys.stderr)
			sys.exit()

		kwargs['combo'] = args.proteins
		kwargs['bins'] = bins
		TM = Crosshits(**kwargs)

	elif args.mode.startswith('cl'):
		kwargs['filename'] = args.clusters
		kwargs['dropclusters'] = args.dropclusters
		kwargs['usegrey'] = args.greyscale
		TM = Clusters(**kwargs)


	if args.show or not args.outfile:
		TM.showTree()

	if args.outfile:
		TM.renderTree(args.outfile, width=args.width, dpi=90)
