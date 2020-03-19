# UROPS

import numpy as np
import pandas as pd
from pandas import Series, DataFrame
import csv

# open the .tsv file
component = DataFrame(csv.reader(open('/Users/suqi/Desktop/Component_6_gene_graphics.tsv'),  delimiter = '\t'))

# Make the first row as the header of the columns.
new_header = component.iloc[0]
component = component[1:]
component.columns = new_header

# delete the empty useless columns
for item in ['FC','SS','Set']:
    del component[item]

# Get the name of the genomes from all the entries
list1 = component['Genome'].tolist()

# Get a non-repetitive list of the names of the genome
list2 = []

for value in list1:
    if value not in list2:
        list2.append(value)

# convert the values in the column "Size (nt)" from strings to numbers
component['Size (nt)'] = pd.to_numeric(component['Size (nt)'])

# define the two possible names of the rSAM protein from the entries
protein1 = 'Radical_SAM'
protein2 = 'SPASM'

# identify the list of genomes with an rSAM protein in the cluster
list3 = []
for organism in list2:
    for index, row in component.iterrows():
        if protein1 or protein2 in row['Function']:
            if organism not in list3:
                list3.append(organism)

# identify list of proteins with an rSAM protein in its cluster and with size less or equal to 150 nucleotides
list4 = []
for organism in list3:
    for index, row in component.iterrows():
        if row['Genome'] == organism:
            if row['Size (nt)'] <= 150 and row['Function'] == 'none':
                l = row.tolist()
                list4.append(l)

# make these potential precursor proteins into a new dataframe
potential = pd.DataFrame(list4, columns = component.columns)

# get the list of genome that these potential proteins belong to
organism_new = []
for index, row in potential.iterrows():
    if row['Genome'] not in organism_new:
        organism_new.append(row['Genome'])

# make two separate lists of genomes whose cluster contain a peptidase or ABC transporter respectively
peptidase = []
transporter = []       
for index, row in component.iterrows():
        if row['Genome'] in organism_new:
            if 'Peptidase' or 'Esterase' or 'Trypsin' in row['Function']:
                if row['Genome'] not in peptidase:
                    peptidase.append(row['Genome'])

    
for index, row in component.iterrows():
        if row['Genome'] in organism_new:
            if 'ABC_tran' in row['Function']:
                if row['Genome'] not in transporter:
                    transporter.append(row['Genome'])

# create another 2 columns in the new dataframe containing the potential proteins indicating the presence/absence of peptidase or transporter in their respective clusters
nan_value = '-'
present = 'Yes'
nonpresent = 'No'
          
potential['Peptidase'] = nan_value
potential['Transporter'] = nan_value
for index, row in potential.iterrows():
    if row['Genome'] in peptidase:
        potential.ix[index,'Peptidase'] = present
    else:
        potential.ix[index,'Peptidase'] = nonpresent

for index, row in potential.iterrows():
    if row['Genome'] in transporter:
        potential.ix[index,'Transporter'] = present
    else:
        potential.ix[index,'Transporter'] = nonpresent

# Get Uniprot ID of these proteins
ID_list = potential['ID'].tolist()

# return the dataframe containing the potential proteins with indication on whether their cluster contain peptidase/transporter
potential
