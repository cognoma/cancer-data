
# coding: utf-8

# # Convert downloaded TCGA datasets into sample × gene matrices

# In[1]:

import os

import pandas


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
# This file contains mutation data (which mutations each sample contains) See the [online documentation](https://genome-cancer.soe.ucsc.edu/proj/site/xena/datapages/?dataset=TCGA.PANCAN.sampleMap/PANCAN_mutation&host=https://tcga.xenahubs.net) for `PANCAN_mutation`.

# In[4]:

path = os.path.join('download', 'PANCAN_mutation.tsv.bz2')
snp_mutation_df = (
    pandas.read_table(path)
    .rename(columns={'sample': 'sample_id'})
)
snp_mutation_df.head(2)


# In[5]:

# The mutation data contains duplicated rows. Disturbing issues like this are
# common in bioinformatics. The following
print('Removing {:.2%} of the mutations, which are duplicates'.format(snp_mutation_df.duplicated().sum() / len(snp_mutation_df)))
snp_mutation_df.drop_duplicates(inplace=True)


# In[6]:

# Number of samples with at least one mutation
snp_mutation_df.sample_id.nunique()


# In[7]:

# Mutations counts by type
snp_mutation_df.effect.value_counts().reset_index()


# ### Convert SNP mutations to gene mutations
# 
# The next cell specifies which mutations to preserve as gene-affecting, which were chosen according to the red & blue [mutation effects in Xena](http://xena.ucsc.edu/how-we-characterize-mutations/).

# In[8]:

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


# In[9]:

# Mutations effects that were observed but nut included
set(snp_mutation_df.effect.unique()) - mutations


# In[10]:

gene_mutation_df = (snp_mutation_df
    .query("effect in @mutations")
    .groupby(['sample_id', 'gene'])
    .apply(len)
    .reset_index()
    .rename(columns={0: 'count'})
)
gene_mutation_df.head(2)


# In[11]:

# Create a sample (rows) by gene (columns) matrix of mutation status
gene_mutation_mat_df = (gene_mutation_df
    .pivot_table(index='sample_id', columns='gene', values='count', fill_value=0)
    .astype(bool).astype(int)
)
gene_mutation_mat_df.shape


# In[12]:

'{:.2%} sample-gene pairs are mutated'.format(
    gene_mutation_mat_df.stack().mean())


# In[13]:

# Top mutated genes
gene_mutation_df.gene.value_counts().reset_index().head(5)


# In[14]:

# Top mutated samples
gene_mutation_df.sample_id.value_counts().reset_index().head(5)


# ## Read gene expression data
# 
# This file contains gene expression data from RNA-Sequencing. See the [online documentation](https://genome-cancer.soe.ucsc.edu/proj/site/xena/datapages/?dataset=TCGA.PANCAN.sampleMap/HiSeqV2&host=https://tcga.xenahubs.net) for `HiSeqV2`.

# In[15]:

# Read the gene × sample dataset
path = os.path.join('download', 'HiSeqV2.tsv.bz2')
expr_df = pandas.read_table(path, index_col=0)

# Process the dataset
expr_df = (expr_df
    # Remove genes containing a `?`
    [~expr_df.index.str.contains('?', regex=False)]
    # Transpose so the data is sample × gene
    .transpose()
    # Sort rows and columns
    .sort_index(axis='rows')
    .sort_index(axis='columns')
)

expr_df.index.rename('sample_id', inplace=True)

expr_df.shape


# In[16]:

# Peak at the data matrix
expr_df.iloc[:5, :5]


# ## Integrate expression and mutation data
# 
# Find samples with both mutation and expression data. We assume that if a sample was not in `PANCAN_mutation`, it was not assayed for mutation. Hence, zero-mutation cancers are excluded even if they have mutation data.

# In[17]:

sample_ids = list(gene_mutation_mat_df.index & expr_df.index)
len(sample_ids)


# In[18]:

# Filter expression (x) and mutation (y) matrices for common samples
x_df = expr_df.loc[sample_ids, :]
y_df = gene_mutation_mat_df.loc[sample_ids, :]


# ### Export matrices to TSVs
# 
# Matrices are saved as sample × gene TSVs. Subsetted matrices are also exported to allow users to quickly explore small portions of the dataset.

# In[19]:

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


# In[20]:

tsv_args = {'sep': '\t', 'float_format': '%.3g'}

for df, name in (x_df, 'expression-matrix'), (y_df, 'mutation-matrix'):

    # Save full dataset
    path = os.path.join('data', name + '.tsv.bz2')
    df.to_csv(path, **tsv_args, compression='bz2')
    
    # Save subsetted datasets
    for sample, nrows, ncols in ('small', 50, 15), ('all-samples', None, 15), ('all-genes', 50, None):
        path = os.path.join('data', 'subset', '{}-{}.tsv'.format(name, sample))
        sample_df(df, nrows=nrows, ncols=ncols).to_csv(path, **tsv_args)

