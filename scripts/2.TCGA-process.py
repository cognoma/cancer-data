
# coding: utf-8

# # Convert downloaded TCGA datasets into sample × gene matrices

# In[1]:

import collections
import os

import pandas


# ## Read gene information

# In[2]:

# Load genes
path = os.path.join('download', 'genes', 'genes.tsv')
gene_df = (pandas.read_table(path)
    .set_index('entrez_gene_id', drop=False)
    [['entrez_gene_id', 'symbol', 'description', 'gene_type']]
)
gene_df.head(2)


# In[3]:

# Load gene updater
path = os.path.join('download', 'genes', 'updater.tsv')
updater_df = pandas.read_table(path)
old_to_new_entrez = dict(zip(updater_df.old_entrez_gene_id, updater_df.new_entrez_gene_id))


# In[4]:

# Load chromosome-symbol to entrez_gene_id mapping
path = os.path.join('download', 'genes', 'chromosome-symbol-mapper.tsv')
chr_sym_map_df = pandas.read_table(path)
chr_sym_map_df.chromosome = 'chr' + chr_sym_map_df.chromosome
chr_sym_map_df.head(2)


# ## Read sample information
# 
# This file contains sample information. See the [online documentation](https://genome-cancer.soe.ucsc.edu/proj/site/xena/datapages/?dataset=TCGA.PANCAN.sampleMap/PANCAN_clinicalMatrix&host=https://tcga.xenahubs.net) for `PANCAN_clinicalMatrix`.
# 
# See [cognoma/cancer-data#14](https://github.com/cognoma/cancer-data/issues/14#issuecomment-238642439 "GitHub Issue: Variable documentation for Xena Browser datasets") for additional variable documentation.

# In[5]:

path = os.path.join('download', 'diseases.tsv')
disease_df = pandas.read_table(path)
disease_df.head(2)


# In[6]:

path = os.path.join('download', 'PANCAN_clinicalMatrix.tsv.bz2')

# Mapping to rename and filter columns
renamer = collections.OrderedDict([
    ('sampleID', 'sample_id'),
    ('_PATIENT', 'patient_id'),
    ('sample_type', 'sample_type'),
    ('_primary_disease', 'disease'),
    ('acronym', 'acronym'),
    ('_primary_site', 'organ_of_origin'),
    ('gender', 'gender'),
    ('age_at_initial_pathologic_diagnosis', 'age_diagnosed'),
    ('_OS_IND', 'dead'),
    ('_OS', 'days_survived'),
    ('_RFS_IND', 'recurred'),
    ('_RFS', 'days_recurrence_free'),
])

clinmat_df = (
    pandas.read_table(path)
    .rename(columns=renamer)
    .merge(disease_df, how='left')
    [list(renamer.values())]
    .set_index('sample_id', drop=False)
    .sort_values('sample_id')
)

# Fix capitalization of gender
clinmat_df.gender = clinmat_df.gender.str.title()

# Save unfiltered dataset to a TSV
path = os.path.join('data', 'complete', 'samples.tsv')
clinmat_df.to_csv(path, sep='\t', float_format='%.0f', index=False)

# Keep only these sample types
# filters duplicate samples per patient
sample_types = {
    'Primary Tumor',
    'Primary Blood Derived Cancer - Peripheral Blood',
}
clinmat_df.query("sample_type in @sample_types", inplace=True)

# Check that no patients are duplicated
assert not clinmat_df.duplicated('patient_id', keep=False).any()

# Check that all diseases in clinmat_df are in disease_df
assert not set(clinmat_df.disease) - set(disease_df.disease)

len(clinmat_df)


# In[7]:

clinmat_df.head(2)


# ## Read mutation data
# 
# This file contains mutation data (which mutations each sample contains) See the [online documentation](https://genome-cancer.soe.ucsc.edu/proj/site/xena/datapages/?dataset=TCGA.PANCAN.sampleMap/PANCAN_mutation&host=https://tcga.xenahubs.net) for `PANCAN_mutation`. Note that duplicate mutation rows, which [occur](https://groups.google.com/d/msg/ucsc-cancer-genomics-browser/eg6nJOFSefw/Z0BM6pU9BAAJ "Message on the Xena Browser Google Group") for samples that were sequenced multiple times, are filtered.

# In[8]:

# Load chromosome-symbol to entrez_gene_id mapping
path = os.path.join('download', 'genes', 'chromosome-symbol-mapper.tsv')
chr_sym_map_df = pandas.read_table(path)
chr_sym_map_df.chromosome = 'chr' + chr_sym_map_df.chromosome
chr_sym_map_df.head(2)


# In[9]:

path = os.path.join('download', 'PANCAN_mutation.tsv.bz2')
snp_mutation_df = (
    pandas.read_table(path)
    .rename(columns={'sample': 'sample_id', 'gene': 'symbol', 'chr': 'chromosome'})
    .merge(chr_sym_map_df)
    .drop(['chromosome', 'symbol'], axis='columns')
    .drop_duplicates()
)
snp_mutation_df.head(2)


# In[10]:

# Number of samples with at least one mutation
snp_mutation_df.sample_id.nunique()


# In[11]:

# Mutations counts by type
snp_mutation_df.effect.value_counts().reset_index()


# ### Convert SNP mutations to gene mutations
# 
# The next cell specifies which mutations to preserve as gene-affecting, which were chosen according to the red & blue [mutation effects in Xena](http://xena.ucsc.edu/how-we-characterize-mutations/).

# In[12]:

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


# In[13]:

# Mutations effects that were observed but nut included
set(snp_mutation_df.effect.unique()) - mutations


# In[14]:

gene_mutation_df = (snp_mutation_df
    .query("effect in @mutations")
    .groupby(['sample_id', 'entrez_gene_id'])
    .apply(len)
    .reset_index()
    .rename(columns={0: 'count'})
)

gene_mutation_df.head(2)


# In[15]:

# Create a sample (rows) by gene (columns) matrix of mutation status
gene_mutation_mat_df = (gene_mutation_df
    .pivot_table(index='sample_id', columns='entrez_gene_id', values='count', fill_value=0)
    .astype(bool).astype(int)
)
gene_mutation_mat_df.shape


# In[16]:

'{:.2%} sample-gene pairs are mutated'.format(
    gene_mutation_mat_df.stack().mean())


# In[17]:

# Save complete mutation matrix
path = os.path.join('data', 'complete', 'mutation-matrix.tsv.bz2')
gene_mutation_mat_df.to_csv(path, sep='\t', compression='bz2')


# ## Read gene expression data
# 
# This file contains gene expression data from RNA-Sequencing. See the [online documentation](https://genome-cancer.soe.ucsc.edu/proj/site/xena/datapages/?dataset=TCGA.PANCAN.sampleMap/HiSeqV2&host=https://tcga.xenahubs.net) for `HiSeqV2`.

# In[18]:

# Read the gene × sample dataset
path = os.path.join('download', 'HiSeqV2.tsv.bz2')
expr_df = pandas.read_table(path, index_col=0)


# In[19]:

# Retrieve symbol to gene mapping for HiSeqV2
path = os.path.join('mapping', 'HiSeqV2-genes', 'HiSeqV2-gene-map.tsv')
gene_map_df = pandas.read_table(path)
symbol_to_entrez = dict(zip(gene_map_df.symbol, gene_map_df.entrez_gene_id))

# Check that there aren't any unmapped symbols
unmapped_symbols = set(expr_df.index) - set(symbol_to_entrez)
unmapped_symbols


# In[20]:

# Process the dataset
expr_df = (expr_df
    # Convert gene symbols to entrez gene ids
    .rename(index=symbol_to_entrez)
    .rename(index=old_to_new_entrez)
    # Average expression for rows with the same entrez_gene_id
    .groupby(level=0).mean()
    # Transpose so the data is sample × gene
    .transpose()
    # Sort rows and columns
    .sort_index(axis='rows')
    .sort_index(axis='columns')
)

expr_df.index.rename('sample_id', inplace=True)

# Filter for valid Entrez Genes
expr_df = expr_df.loc[:, expr_df.columns.isin(gene_df.entrez_gene_id)]

expr_df.shape


# In[21]:

# Number of patients represented in the expression dataset
clinmat_df.query("sample_id in @expr_df.index").patient_id.nunique()


# In[22]:

# Peak at the data matrix
expr_df.iloc[:5, :5]


# In[23]:

# Save complete expression matrix
path = os.path.join('data', 'complete', 'expression-matrix.tsv.bz2')
expr_df.to_csv(path, sep='\t', float_format='%.3g', compression='bz2')


# ## Integrate expression and mutation data
# 
# Find samples with both mutation and expression data. We assume that if a sample was not in `PANCAN_mutation`, it was not assayed for mutation. Hence, zero-mutation cancers are excluded even if they have mutation data.

# In[24]:

sample_ids = sorted(clinmat_df.index & gene_mutation_mat_df.index & expr_df.index)
len(sample_ids)


# In[25]:

# Filter expression (x) and mutation (y) matrices for common samples
sample_df = clinmat_df.loc[sample_ids, :]
x_df = expr_df.loc[sample_ids, :]
y_df = gene_mutation_mat_df.loc[sample_ids, :]


# In[26]:

# Add a columnn for mutations per sample
sample_df['n_mutations'] = y_df.sum(axis='columns')


# In[27]:

x_gene_df = (
    gene_df.loc[x_df.columns, :]
    .assign(mean_expression=x_df.mean(axis='rows'))
    .assign(stdev_expression=x_df.std(axis='rows'))
)
path = os.path.join('data', 'expression-genes.tsv')
x_gene_df.to_csv(path, sep='\t', index=False, float_format='%.4g')
x_gene_df.head(2)


# In[28]:

y_gene_df = (
    gene_df.loc[y_df.columns, :]
    .assign(n_mutations=y_df.sum(axis='rows'))
    .assign(mutation_freq=y_df.mean(axis='rows'))
)
path = os.path.join('data', 'mutation-genes.tsv')
y_gene_df.to_csv(path, sep='\t', index=False, float_format='%.4g')
y_gene_df.head(2)


# ### Export matrices to TSVs
# 
# Matrices are saved as sample × gene TSVs. Subsetted matrices are also exported to allow users to quickly explore small portions of the dataset.

# In[29]:

path = os.path.join('data', 'samples.tsv')
sample_df.to_csv(path, sep='\t', float_format='%.0f', index=False)


# In[30]:

def subset_df(df, nrows=None, ncols=None, row_seed=0, col_seed=0):
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


# In[31]:

tsv_args = {'sep': '\t', 'float_format': '%.3g'}

for df, name in (x_df, 'expression-matrix'), (y_df, 'mutation-matrix'):

    # Save full dataset
    path = os.path.join('data', name + '.tsv.bz2')
    df.to_csv(path, **tsv_args, compression='bz2')
    
    # Save subsetted datasets
    for sample, nrows, ncols in ('small', 50, 15), ('all-samples', None, 15), ('all-genes', 50, None):
        path = os.path.join('data', 'subset', '{}-{}.tsv'.format(name, sample))
        subset_df(df, nrows=nrows, ncols=ncols).to_csv(path, **tsv_args)

