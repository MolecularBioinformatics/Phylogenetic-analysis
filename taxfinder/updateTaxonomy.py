#!/usr/bin/env python3

import sys

if sys.argv[1] == 'acc2taxid':
	with open('unsorted_acc2taxid', 'w') as out:
		print('py: reading pdb_chain_taxonomy.tsv and writing unsorted_acc2taxid...')
		with open('pdb_chain_taxonomy.tsv', 'r') as f:
			next(f)
			next(f)
			pdb = None
			taxid = None
			for line in f:
				lline = line.split()
				if pdb == lline[0]:
					continue
				if pdb is not None:
					out.write('{:<12}{:<7}\n'.format(pdb.upper(), taxid))
				pdb = lline[0]
				taxid = lline[2]
		print('py: reading prot.accession2taxid and writing unsorted_acc2taxid...')
		with open('prot.accession2taxid', 'r') as f:
			next(f)
			for line in f:
				lline = line.split()
				# Writes only the accession and the taxid to the new file.
				# Accession and taxid are right-padded with spaces:
				# acc12345____1234___{NEWLINE}
				# ^--- 12 ---^^- 7 -^
				# This allows later for an easier binary search in the file.
				# A field separator (like tab) is not necessary as the column width is constant
				out.write('{:<12}{:<7}\n'.format(lline[0], lline[2]))

if sys.argv[1] == 'taxinfo':
	data = {}	# data[taxid] = [depth, parent, rank, name]

	print('py: reading names.dmp...')

	with open('names.dmp', 'r') as f:
		for line in f:
			lline = line.split('\t|\t')
			# taxid, name, uniqueName, nameClass
			if lline[3] == 'scientific name\t|\n':
				data[lline[0]] = ['@', '@', '@', lline[1]]

	print('py: reading nodes.dmp...')

	with open('nodes.dmp', 'r') as f:
		for line in f:
			lline = line.split('\t|\t')
			# taxid, parentTaxid, rank, *others
			data[lline[0]][1] = lline[1]
			data[lline[0]][2] = lline[2]

	print('py: writing taxinfo...')

	with open('taxinfo', 'w') as out:
		for taxid in sorted(data.keys()):
			# TaxID, Level, Parent, Rank, Name
			level = 0
			tid = taxid
			while tid != '1':
				tid = data[tid][1]
				level += 1
			data[taxid][0] = str(level)
			out.write('{}\t{}\n'.format(taxid, '\t'.join(data[taxid])))
