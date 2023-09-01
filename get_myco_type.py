"""
Looks up mycorrhizal type for a given species (ID'd using 
genus and epithet), and falls back onto the family if no
data is found for this species in our databases

Databases (in order in which we look):
* Fungalroot
* MycoDB

Used as a command line program

Needs:
* Pandas to be installed (install with "pip install pandas" or "pip3 install pandas")
"""

import pandas as pd

# load MycoDB
mycodb_csv = pd.read_csv('mycodb/MycoDB_version4.csv')

def lookup_mycodb(family, out=False):
    # PlantSpecies
    rows = mycodb_csv['PlantFamily'] == family.lower()
    
    type_count = dict()
    
    for index, row in mycodb_csv[rows].iterrows():
        mtype = row['MYCORRHIZAETYPE']
        
        if mtype not in type_count:
            type_count[mtype] = 0
        
        type_count[mtype] += 1
    
    if out and len(type_count) > 1:
        print('warning: number of mycorrhizae types present at this plant species is', len(type_count), 'species:', type_count)
    
    mtype = None
    best_count = 0
    
    for key in type_count:
        if type_count[key] > best_count:
            best_count = type_count[key]
            mtype = key
    
    if out and mtype is None:
        print('Warning: no data found for this plant family:', family.lower())
    
    return mtype

# load fungalroot
fungalrootdb = (pd.read_csv('fungalroot/occurrence.csv'), pd.read_csv('fungalroot/measurements.csv'), pd.read_csv('fungalroot/appendix.csv'))

def lookup_fungalroot(family, species):
    parts = species.split()
    genus = parts[0]
    epithet = parts[1]
    
    rows = fungalrootdb[0][fungalrootdb[0]['genus'] == genus]
    
    id_ = None
    
    # assume there is only one occurrence of the plant in occurrences.txt
    for index, row in rows.iterrows():
        if row.specificEpithet == epithet:
            id_ = row.identifier
            break
    
    rows = fungalrootdb[1][fungalrootdb[1]['Core ID'] == id_]
    
    for index, row in rows.iterrows():
        if row.measurementType == 'Mycorrhiza type':
            return row.measurementValue
            break
    
    stuff = fungalrootdb[2].loc[fungalrootdb[2]['Genus'] == genus]['Mycorrhizal type']
    
    if len(stuff) != 1:
        return None
    
    return stuff.item()

# function that checks fungalroot, then mycodb
def get_myco_type(family, species):
    x = lookup_fungalroot(family, species)
    
    if x == 'non-mycorrhizal':
        x = 'NM'
    
    if x is not None:
        return x
    
    return lookup_mycodb(family)

print('Welcome to the mycorrhizal type lookup tool! This tool will search Fungalroot, then MycoDB, to see ' + 
	'if we can find a mycorrhizal type associated with a given species. If a species is not in our databases, ' +
	'we will attempt to find the mycorrhizal type associated with its family instead. ')
print('Example: Liliaceae Trillium cernuum L. => AM')
print('This is what it might look like:')
print('What is this family (if known)? If not known press Enter for blank: Liliaceae')
print('What is the species: Trillium cernuum L.')
print('AM')
print()

print('----------------------------')

while True:
	family = input('What is this family (if known)? If not known press Enter for blank: ')
	species = input('What is the species: ')

	print(get_myco_type(family, species))