
# coding: utf-8

# ### Convert mutation gene symbol labels to Entrez IDs  
# Goal: Relabel the mutation data frame with Entrez IDs instead of gene names, by mapping a combination of chromosome and gene symbol to Entrez ID. The NCBI file downloaded and read in the next cell contains the Entrez ID - gene symbol pairs we will use to do so.

# In[1]:

import os
import pandas
from urllib.request import urlretrieve
import numpy as np


# ### Creating Mappings

# In[2]:

base_url = 'ftp://ftp.ncbi.nih.gov/gene/DATA/GENE_INFO/Mammalia/'
map_name = 'Homo_sapiens.gene_info.gz'
map_url = base_url + map_name
path = os.path.join('download', map_name)
urlretrieve(map_url, path)

path = os.path.join('download', map_name)

map_df = (
    pandas.read_table(path)
    .rename(columns={'#tax_id' : 'tax_id', 'GeneID': 'entrez_id', 'chromosome':'chr', 'Symbol':'symbol'})
)

#filter to include only tax_id == 9606 (human) and columns of interest
map_df = map_df.ix[map_df['tax_id'] == 9606, ['entrez_id','chr', 'symbol', 'Synonyms']]

map_df.head(2)


# In[3]:

map_primary = map_df[['entrez_id','chr', 'symbol']]
map_primary.head(2)


# In[4]:

alternates = map_df['Synonyms'].str.split('|').apply(pandas.Series, 1).stack()
alternates.name = 'alt_symbol' # needs a name to join
alternates.index = alternates.index.droplevel(-1)
map_alternates = map_df.join(alternates)


# In[5]:

#format as entrez_id, chr, alternate symbol
map_alternates = map_alternates.drop_duplicates(subset=['chr','alt_symbol'], keep=False)[
    ['entrez_id', 'chr', 'alt_symbol']].rename(columns={'alt_symbol': 'symbol'})


# In[6]:

#use keep = first to give primacy to the primary symbol convention in the case of ambiguous mappings
maps_combined = map_primary.append(map_alternates).drop_duplicates(subset=['chr','symbol'], keep='first')

#add chr string to facilitate integration with mutation dataset
maps_combined['chr'] = 'chr' + maps_combined['chr'].astype(str)
maps_combined.head(2)


# In[7]:

# Check that chr/symbol are all unique 
assert not (maps_combined.chr + '|' + maps_combined.symbol).duplicated().any()


# ### Check and see which mutations in the dataset fail to map

# In[8]:

path = "../../download/PANCAN_mutation.tsv.bz2"
mutation_df = pandas.read_table(path, index_col=0)

mutation_df.head(2)


# In[9]:

failed_mappings = (set(mutation_df.chr + '|' + mutation_df.gene) - 
                   set(maps_combined.chr + '|' + maps_combined.symbol))
    
'{0} of {1} mutations failed to map based on chromosome and either primary or alternate gene symbol. ({2:.2%} of mutations.)'.format(
    len(failed_mappings),
    len(mutation_df.chr),
    len(failed_mappings)/len(mutation_df.chr))


# Some (~300) of these failed mappings are attributable to non-standard chromosomes designations. 

# In[10]:

#count the mutations observed on chromosomes that failed to map for any observations
failed_chr_mappings = (set(mutation_df.chr ) - set(maps_combined.chr))

pandas.merge(mutation_df, pandas.DataFrame(list(failed_chr_mappings)), 
             left_on=['chr'], right_on=[0] , how='inner').chr.value_counts()


# In[11]:

#remove the alternate gene symbols 'NaN' corresponding to entrez id 11280  and 'NA' correspongind to 7504
maps_combined = maps_combined[(maps_combined.symbol != 'NaN') & (maps_combined.symbol != 'NA')]

maps_combined.shape


# ### Export Mappings

# In[12]:

maps_combined.to_csv('PANCAN-mutation-gene-map.tsv', index=False, sep='\t')


# In[ ]:



