#UROPS
import numpy as np
import pandas as pd
from pandas import Series, DataFrame
import csv
import sys

#open the .tsv file. Run in terminal. When in the working directory, input in the terminal: python #program_file_name #component_number
filename = sys.argv[1]

component = DataFrame(csv.reader(open('/Users/suqi/Desktop/Component_%s_gene_graphics.tsv' % (filename)),  delimiter = '\t')) # change the directory to where the tsv file is

#Make the first row as the header of the columns. And delete the empty useless columns
new_header = component.iloc[0]
component = component[1:]
component.columns = new_header

for item in ['FC','SS','Set']:
    del component[item]

#Get the name of the genomes from all the entries and convert to a non-repetitive list
list1 = component['Genome'].tolist()

list2 = []

for value in list1:
    if value not in list2:
        list2.append(value)

#convert the values in the column "Size (nt)" from strings to numbers
component['Size (nt)'] = pd.to_numeric(component['Size (nt)'])

#define the two possible names of the rSAM protein from the entries and identify the list of genomes with an rSAM protein in the cluster
protein1 = 'Radical_SAM'
protein2 = 'SPASM'

list3 = []
for organism in list2:
    for index, row in component.iterrows():
        if protein1 or protein2 in row['Function']:
            if organism not in list3:
                list3.append(organism)

#identify list of proteins with an rSAM protein in its cluster and with size less or equal to 150 nucleotides
list4 = []
for organism in list3:
    for index, row in component.iterrows():
        if row['Genome'] == organism:
            if row['Size (nt)'] <= 150 and row['Function'] == 'none':
                l = row.tolist()
                list4.append(l)

#make these potential precursor proteins into a new dataframe and get the list of genome that these potential proteins belong to
potential = pd.DataFrame(list4, columns = component.columns)

organism_new = []
for index, row in potential.iterrows():
    if row['Genome'] not in organism_new:
        organism_new.append(row['Genome'])

#make two separate lists of genomes whose cluster contain a peptidase or ABC transporter respectively
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

#create another 2 columns in the new dataframe containing the potential proteins indicating the presence/absence of peptidase or transporter in their respective clusters
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

#Get Uniprot ID of these proteins and return the dataframe containing the potential proteins with indication on whether their cluster contain peptidase/transporter
ID_list = potential['ID'].tolist()

potential

#put all the id into a string separated by space.
ID = ""

for item in ID_list:
    ID = ID + " " + item
    
print(ID)

from Bio import ExPASy
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_protein

#programmatically retrieve the data from Uniprot into a list (seqlist) of sequence records. When encounter an error (eg. protein record deleted from Uniprot), skip that record while printing the error so it will ot crash
n = 1
seqlist = []
for uniprot_id in ID.split():
    try:
        with ExPASy.get_sprot_raw(uniprot_id) as handle:
            seq_record = SeqIO.read(handle, "swiss")
            print("seq ", n)
            print(seq_record.id)
            print(seq_record.name)
            print(seq_record.description)
            print(repr(seq_record.seq))
            print("Length %i" % len(seq_record))
            n = n + 1
            seqlist.append(seq_record)
    except Exception as e:
        print("Error: ", e , uniprot_id)
        
#write the sequence record of these proteins into a fasta file in the working directory that can be referred to if needed
SeqIO.write(seqlist, "component_%s.fasta" % (filename), "fasta")

#Identify the protein seqeunces with the known pattern for the precursor: WSW or WNW in the 10 amino acid from C terminal.This list: precursor_1 contain the most probable precursors since they have the known pattern.
precursor_1 = []

for record in seqlist:
    last_ten_aa = str(record.seq[-10:])
    if "WSW" in last_ten_aa:
        precursor_1.append(record.id)
    elif "WNW" in last_ten_aa:
        precursor_1.append(record.id)
        
        
#Identify the sequences without the known pattern but contain more than 1 aromatic amino acid in the 10 aa from C terminal. Put their ID in the prob_precursor
prob_precursor = []

for record in seqlist:
    if record not in precursor_1:
        last_ten_aa = str(record.seq[-10:])
        aro_freq = 0
        for item in ["F", "W", "Y", "H"]:
            item_freq = last_ten_aa.count(item)
            aro_freq = aro_freq + item_freq
        if aro_freq > 1:
            prob_precursor.append(record.id)


#for those in the probable precursors, check in the dataframe formed above for whether they have the two characteristic proteins (transporter and protease) in its cluster
prob_precursor_1 = [] #have both
prob_precursor_2 = [] #have one
prob_precursor_3 = [] #have none

for index, row in potential.iterrows():
    if row['ID'] in prob_precursor:
        if row['Peptidase'] == "Yes":
            if row['Transporter'] == "Yes":
                prob_precursor_1.append(row['ID'])
            else:
                prob_precursor_2.append(row['ID'])
        else:
            if row['Transporter'] == "Yes":
                prob_precursor_2.append(row['ID'])
            else:
                prob_precursor_3.append(row['ID'])
                
print("Both: ", prob_precursor_1)
print("One: ", prob_precursor_2)
print("None: ", prob_precursor_3)

#Get other patterns, train the code to recognize 
