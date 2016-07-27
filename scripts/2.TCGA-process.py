
# coding: utf-8

# # Convert downloaded TCGA datasets into sample × gene matrices

# In[1]:

import os
import pandas
from urllib.request import urlretrieve
import numpy as np


# ## Read sample information
# 
# This file contains sample information. See the [online documentation](https://genome-cancer.soe.ucsc.edu/proj/site/xena/datapages/?dataset=TCGA.PANCAN.sampleMap/PANCAN_clinicalMatrix&host=https://tcga.xenahubs.net) for `PANCAN_clinicalMatrix`.

# In[2]:

path = os.path.join('download', 'PANCAN_clinicalMatrix.tsv.bz2')
clinmat_df = (
    pandas.read_table(path)
    .rename(columns={'sampleID': 'sample_id'})
)
# Check that no sample_ids are duplicated
assert not clinmat_df.sample_id.duplicated().any()
clinmat_df.shape


# In[3]:

# Types of samples
clinmat_df.sample_type.value_counts()


# ## Read mutation data
# 
# This file contains mutation data (which mutations each sample contains) See the [online documentation](https://genome-cancer.soe.ucsc.edu/proj/site/xena/datapages/?dataset=TCGA.PANCAN.sampleMap/PANCAN_mutation&host=https://tcga.xenahubs.net) for `PANCAN_mutation`. Note that duplicate mutation rows, which [occur](https://groups.google.com/d/msg/ucsc-cancer-genomics-browser/eg6nJOFSefw/Z0BM6pU9BAAJ "Message on the Xena Browser Google Group") for samples that were sequenced multiple times, are filtered.

# In[4]:

path = os.path.join('download', 'PANCAN_mutation.tsv.bz2')
snp_mutation_df = (
    pandas.read_table(path)
    .rename(columns={'sample': 'sample_id'})
    .drop_duplicates()
)
snp_mutation_df.head(2)


# In[5]:

# Number of samples with at least one mutation
snp_mutation_df.sample_id.nunique()


# In[6]:

# Mutations counts by type
snp_mutation_df.effect.value_counts().reset_index()


# ### Convert mutation gene symbol labels to Entrez IDs  

# Goal: Relabel the mutation data frame with Entrez IDs instead of gene names, by mapping a combination of chromosome and gene symbol to Entrez ID. The NCBI file downloaded and read in the next cell contains the Entrez ID - gene symbol pairs we will use to do so.

# In[7]:

base_url = 'ftp://ftp.ncbi.nih.gov/gene/DATA/GENE_INFO/Mammalia/'
map_name = 'Homo_sapiens.gene_info.gz'
map_url = base_url + map_name
path = os.path.join('download', map_name)
urlretrieve(map_url, path)

path = os.path.join('download', 'Homo_sapiens.gene_info.gz')
map_names = ['tax_id', 'GeneID', 'Symbol', 'LocusTag', 'Synonyms', 'dbXrefs',
             'chromosome', 'map_location', 'description', 'type_of_gene',
             'Symbol_from_nomenclature_authority', 'Full_name_from_nomenclature_authority',
             'Nomenclature_status', 'Other_designations', 'Modification_date']

#skip the first row, which contains a description of the headers
map_df = (
    pandas.read_table(path,  names = map_names, skiprows = 1)
    .rename(columns={'GeneID': 'entrez_id', 'chromosome':'chr'})
)

map_df.head(2)


# First, we will map base on unambiguous combinations of chromosome and the gene symbol of record.

# In[8]:

#create chr-gene designation and check for duplicates
snp_mutation_df['chr_gene'] = (snp_mutation_df['chr'] + '-' + snp_mutation_df['gene'])
map_df['chr_gene'] = 'chr' + map_df['chr'] + '-' + map_df['Symbol']

'{0} of the {1} possible chr-gene combinations in map_df are unique.'.format(
    map_df['chr_gene'].nunique(),
    len(map_df['chr_gene'])
    )


# In[9]:

#remove all duplicated chr-gene combinations to avoid ambiguous mapping
map_df_nodups = map_df.drop_duplicates(subset='chr_gene', keep=False)

'{0} of the {1} chr-gene combinations in map_df_nodups are unique.'.format(
    map_df_nodups['chr_gene'].nunique(),
    len(map_df_nodups['chr_gene'])
    )


# In[10]:

#perform intial merge
map_entrez = pandas.merge(snp_mutation_df, map_df_nodups, on = 'chr_gene', how = 'inner', suffixes = ['','_y'])

#assess merge completeness
all_maps = pandas.merge(snp_mutation_df, map_df_nodups, on = 'chr_gene', how = 'left', suffixes=['', '_y'])
fail_index = all_maps.index[all_maps['entrez_id'].apply(np.isnan)]
fail_maps = all_maps.iloc[fail_index]

#drop the unsuccessfully merged columns from the unmapped observations
fail_maps = fail_maps[snp_mutation_df.columns.values]

'Attempted to map {0} total observations. Mapped {1} observations unambiguously. Failed to map {2} observations. {3:.2%} of the mutations failed to map based on chr + gene symbol'.format(
snp_mutation_df.shape[0],
map_entrez.shape[0],
fail_maps.shape[0],
fail_maps.shape[0]/(map_entrez.shape[0]+fail_maps.shape[0] )  
)


# Attempt to use alternate gene symbols to re-map observations that intially failed to map.

# In[11]:

#extract alternate symbol info, create row for each alternate
alternates = map_df['Synonyms'].str.split('|').apply(pandas.Series, 1).stack()
alternates.name = 'alt_symbol' # needs a name to join
alternates.index = alternates.index.droplevel(-1)
map_alternates = map_df.join(alternates)

map_alternates['chr_gene'] = ('chr' + map_alternates['chr'] + '-' + map_alternates['alt_symbol'])

#remove all duplicated chr-gene combinations to avoid ambiguous mapping
map_alternates_nodups = map_alternates.drop_duplicates(subset='chr_gene', keep=False)

#merge and completeness based on alternate symbols
map_entrez_alts = pandas.merge(fail_maps, map_alternates_nodups, on = 'chr_gene', how = 'inner', suffixes = ['', '_y'])

'An additional {0} observations mapped based on chr + alternate gene symbols. This represents {1:.2%} of observations that initially failed to map. {2:.2%} of total observations remain unmapped. '.format(
    map_entrez_alts.shape[0],
    map_entrez_alts.shape[0]/fail_maps.shape[0],
    1 - (map_entrez_alts.shape[0] + map_entrez.shape[0])/snp_mutation_df.shape[0],
)


# Combine the Entrez-labelled dataframes for observations mapped based on standard or alternate gene symbols. Keep only entrez_id and columns in the original mutation df.

# In[12]:

cols = np.append(snp_mutation_df.columns.values, 'entrez_id')
snp_mutation_df_mapped = pandas.concat([map_entrez[cols], map_entrez_alts[cols]])
snp_mutation_df_mapped.head(2)


# ### Convert SNP mutations to gene mutations
# 
# The next cell specifies which mutations to preserve as gene-affecting, which were chosen according to the red & blue [mutation effects in Xena](http://xena.ucsc.edu/how-we-characterize-mutations/).

# In[13]:

mutations = {
    'Frame_Shift_Del',
    'Frame_Shift_Ins',
    'In_Frame_Del',
    'In_Frame_Ins',
    'Missense_Mutation',
    'Nonsense_Mutation',
    'Nonstop_Mutation',
    'RNA',
    'Splice_Site',
    'Translation_Start_Site',
}


# In[14]:

# Mutations effects that were observed but not included
set(snp_mutation_df_mapped.effect.unique()) - mutations


# In[15]:

gene_mutation_df = (snp_mutation_df_mapped
    .query("effect in @mutations")
    .groupby(['sample_id', 'entrez_id'])
    .apply(len)
    .reset_index()
    .rename(columns={0: 'count'})
)
gene_mutation_df.head(2)


# In[16]:

# Create a sample (rows) by gene (columns) matrix of mutation status
gene_mutation_mat_df = (gene_mutation_df
    .pivot_table(index='sample_id', columns='entrez_id', values='count', fill_value=0)
    .astype(bool).astype(int)
)
gene_mutation_mat_df.shape


# In[17]:

'{:.2%} sample-gene pairs are mutated'.format(
    gene_mutation_mat_df.stack().mean())


# In[18]:

# Top mutated genes, with table relabelled for clarity
gene_mutation_df.entrez_id.value_counts().reset_index().rename(columns={'index':'entrez_id', 'entrez_id':'count'}).head(5)


# The top mutated gene (Entrez ID: 7157) is TP53.

# In[19]:

# Top mutated samples
gene_mutation_df.sample_id.value_counts().reset_index().head(5)


# ## Read gene expression data
# 
# This file contains gene expression data from RNA-Sequencing. See the [online documentation](https://genome-cancer.soe.ucsc.edu/proj/site/xena/datapages/?dataset=TCGA.PANCAN.sampleMap/HiSeqV2&host=https://tcga.xenahubs.net) for `HiSeqV2`.

# In[20]:

# Read the gene × sample dataset
path = os.path.join('download', 'HiSeqV2.tsv.bz2')
expr_df = pandas.read_table(path, index_col=0)


# In[21]:

# Retrieve symbol to gene mapping for HiSeqV2
path = os.path.join('mapping', 'HiSeqV2-genes', 'HiSeqV2-gene-map.tsv')
gene_map_df = pandas.read_table(path)
symbol_to_entrez = dict(zip(gene_map_df.symbol, gene_map_df.entrez_gene_id))

# Check that there aren't any unmapped symbols
unmapped_symbols = set(expr_df.index) - set(symbol_to_entrez)
unmapped_symbols


# In[22]:

# Process the dataset
expr_df = (expr_df
    # Convert gene symbols to entrez gene ids
    .rename(index=symbol_to_entrez)
    # Transpose so the data is sample × gene
    .transpose()
    # Sort rows and columns
    .sort_index(axis='rows')
    .sort_index(axis='columns')
)

expr_df.index.rename('sample_id', inplace=True)

expr_df.shape


# In[23]:

# Peak at the data matrix
expr_df.iloc[:5, :5]


# ## Integrate expression and mutation data
# 
# Find samples with both mutation and expression data. We assume that if a sample was not in `PANCAN_mutation`, it was not assayed for mutation. Hence, zero-mutation cancers are excluded even if they have mutation data.

# In[24]:

sample_ids = list(gene_mutation_mat_df.index & expr_df.index)
len(sample_ids)


# In[25]:

# Filter expression (x) and mutation (y) matrices for common samples
x_df = expr_df.loc[sample_ids, :]
y_df = gene_mutation_mat_df.loc[sample_ids, :]


# ### Export matrices to TSVs
# 
# Matrices are saved as sample × gene TSVs. Subsetted matrices are also exported to allow users to quickly explore small portions of the dataset.

# In[26]:

def sample_df(df, nrows=None, ncols=None, row_seed=0, col_seed=0):
    """Randomly subset a dataframe, preserving row and column order."""
    if nrows is None:
        nrows = len(df)
    if ncols is None:
        ncols = len(df.columns)
    return (df
        .sample(n=nrows, random_state=row_seed, axis='rows')
        .sample(n=ncols, random_state=col_seed, axis='columns')
        .sort_index(axis='rows')
        .sort_index(axis='columns')
    )


# In[27]:

tsv_args = {'sep': '\t', 'float_format': '%.3g'}

for df, name in (x_df, 'expression-matrix'), (y_df, 'mutation-matrix'):

    # Save full dataset
    path = os.path.join('data', name + '.tsv.bz2')
    df.to_csv(path, **tsv_args, compression='bz2')
    
    # Save subsetted datasets
    for sample, nrows, ncols in ('small', 50, 15), ('all-samples', None, 15), ('all-genes', 50, None):
        path = os.path.join('data', 'subset', '{}-{}.tsv'.format(name, sample))
        sample_df(df, nrows=nrows, ncols=ncols).to_csv(path, **tsv_args)


# In[ ]:



