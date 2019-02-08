#!/usr/bin/env bash

echo "Downloading accession to taxid file..."
curl -R --retry 1 -o prot.accession2taxid.gz ftp://ftp.ncbi.nih.gov/pub/taxonomy/accession2taxid/prot.accession2taxid.gz

echo "Downloading PDB to taxid file..."
curl -R --retry 1 -o pdb_chain_taxonomy.tsv.gz ftp://ftp.ebi.ac.uk/pub/databases/msd/sifts/flatfiles/tsv/pdb_chain_taxonomy.tsv.gz

echo "Unzipping accession to taxid file..."
gunzip prot.accession2taxid.gz

echo "Unzipping PDB to taxid file..."
gunzip pdb_chain_taxonomy.tsv.gz

echo "Running Python to create the acc2taxid database..."
python3 updateTaxonomy.py acc2taxid

echo "Sorting accession to taxid file..."
sort -s -k 1,1 --parallel=6 -o acc2taxid unsorted_acc2taxid

echo "Counting lines of the acc2taxid file..."
wc -l < acc2taxid > numLines

echo "Cleaning up the acc2taxid part..."
rm prot.accession2taxid
rm pdb_chain_taxonomy.tsv
rm unsorted_acc2taxid

echo "Downloading lineage file..."
curl -R --retry 1 -o taxdump.tar.gz ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz

echo "Unzipping lineage file..."
tar xzf taxdump.tar.gz
rm taxdump.tar.gz

echo "Running Python to create the taxonomy database..."
python3 updateTaxonomy.py taxinfo

echo "Cleaning up..."
rm *.dmp
rm gc.prt
rm readme.txt

echo "Databases are created. Now it's time for some holidays..."
