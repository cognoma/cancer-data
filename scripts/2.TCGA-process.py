
# coding: utf-8

# # Convert downloaded TCGA datasets into sample × gene matrices
# 
# This notebook is updated to include the data from the [TCGA PanCanAtlas April 2018 updates](http://www.cell.com/pb-assets/consortium/pancanceratlas/pancan/index.html).

# In[1]:


import collections
import os

import pandas


# ## Read gene information

# In[2]:


# Load genes
path = os.path.join('download', 'genes', 'genes.tsv')
gene_df = (pandas.read_table(path, dtype='str')
    .set_index('entrez_gene_id', drop=False)
    [['entrez_gene_id', 'symbol', 'description', 'gene_type']]
)
gene_df.head(2)


# In[3]:


# Load gene updater
path = os.path.join('download', 'genes', 'updater.tsv')
updater_df = pandas.read_table(path)
old_to_new_entrez = dict(zip(updater_df.old_entrez_gene_id,
                             updater_df.new_entrez_gene_id))


# In[4]:


# Load chromosome-symbol to entrez_gene_id mapping
path = os.path.join('download', 'genes', 'chromosome-symbol-mapper.tsv')
chr_sym_map_df = pandas.read_table(path)
chr_sym_map_df.chromosome = 'chr' + chr_sym_map_df.chromosome
chr_sym_map_df.head(2)


# ## Read sample information
# 
# This file contains sample information. See the [online documentation](https://xenabrowser.net/datapages/?dataset=EB%2B%2BAdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.xena&host=https%3A%2F%2Fpancanatlas.xenahubs.net) for more information. For more details on curation refer to [Liu et al. 2018](https://doi.org/10.1016/j.cell.2018.02.052 "An Integrated TCGA Pan-Cancer Clinical Data Resource to Drive High-Quality Survival Outcome Analytics")

# In[5]:


path = os.path.join('data', 'diseases.tsv')
disease_df = pandas.read_table(path)
disease_df.head(2)


# In[6]:


# Data from https://gdc.cancer.gov/resources-tcga-users/tcga-code-tables/sample-type-codes
path = os.path.join('mapping', 'tcga_sampletype_codes.csv')
sampletype_codes_df = pandas.read_csv(path, dtype='str')
sampletype_codes_dict = dict(zip(sampletype_codes_df.Code, sampletype_codes_df.Definition))
sampletype_codes_df.head(2)


# In[7]:


path = os.path.join('download', 'Survival_SupplementalTable_S1_20171025_xena_sp.tsv.gz')

# Mapping to rename and filter columns
renamer = collections.OrderedDict([
    ('sample', 'sample_id'),
    ('_PATIENT', 'patient_id'),
    ('cancer type abbreviation', 'acronym'),
    ('__placeholder__', 'disease'),
    ('age_at_initial_pathologic_diagnosis', 'age_diagnosed'),
    ('gender', 'gender'),
    ('race', 'race'),
    ('ajcc_pathologic_tumor_stage', 'ajcc_stage'),
    ('clinical_stage', 'clinical_stage'),
    ('histological_type', 'histological_type'),
    ('histological_grade', 'histological_grade'),
    ('initial_pathologic_dx_year', 'initial_pathologic_dx_year'),
    ('menopause_status', 'menopause_status'),
    ('birth_days_to', 'birth_days_to'),
    ('vital_status', 'vital_status'),
    ('tumor_status', 'tumor_status'),
    ('last_contact_days_to', 'last_contact_days_to'),
    ('death_days_to', 'death_days_to'),
    ('cause_of_death', 'cause_of_death'),
    ('new_tumor_event_type', 'new_tumor_event_type'),
    ('new_tumor_event_site', 'new_tumor_event_site'),
    ('new_tumor_event_site_other', 'new_tumor_event_site_other'),
    ('new_tumor_event_dx_days_to', 'days_recurrence_free'),
    ('treatment_outcome_first_course', 'treatment_outcome_first_course'),
    ('margin_status', 'margin_status'),
    ('residual_tumor', 'residual_tumor'),
    ('_EVENT', 'event_status'),
    ('_TIME_TO_EVENT', 'event_days'),
    ('OS', 'dead'),
    ('OS.time', 'days_survived'),
    ('DSS', 'disease_specific_survival_status'),
    ('DSS.time', 'disease_specific_survival_days'),
    ('DFI', 'disease_free_interval_status'),
    ('DFI.time', 'disease_free_interval_days'),
    ('PFI', 'progression_free_interval_status'),
    ('PFI.time', 'progression_free_interval_days'),
    ('Redaction', 'Redaction')
])

clinmat_df = (
    pandas.read_table(path)
    .rename(columns=renamer)
    .merge(disease_df, how='left')
    [list(renamer.values())]
    .set_index('sample_id', drop=False)
    .sort_values('sample_id')
)

# Fix capitalization of gender and race
clinmat_df.gender = clinmat_df.gender.str.title()
clinmat_df.race = clinmat_df.race.str.title()

# Extract sample-type with the code dictionary
clinmat_df = clinmat_df.assign(sample_type = clinmat_df.sample_id.str[-2:])
clinmat_df.sample_type = clinmat_df.sample_type.replace(sampletype_codes_dict)


# In[8]:


clinmat_df.head(2)


# In[9]:


# What is the distribution of cancer-types
clinmat_df.sample_type.value_counts()


# In[10]:


# Save unfiltered dataset to a TSV
path = os.path.join('data', 'complete', 'samples.tsv')
clinmat_df.to_csv(path, sep='\t', float_format='%.0f', index=False)


# In[11]:


# Remove "Redacted" tumors
# These patients either withdrew consent or had genome data mismatch errors
clinmat_df = clinmat_df.query('Redaction != "Redacted"').drop("Redaction", axis=1)

# Keep only these sample types
# filters duplicate samples per patient
sample_types = {
    'Primary Solid Tumor',
    'Primary Blood Derived Cancer - Peripheral Blood',
}
clinmat_df.query("sample_type in @sample_types", inplace=True)


# In[12]:


# Check that no samples are duplicated
assert not clinmat_df.duplicated('sample_id', keep=False).any()

# Check that all diseases in clinmat_df are in disease_df
assert not set(clinmat_df.acronym) - set(disease_df.acronym)

len(clinmat_df)


# In[13]:


clinmat_df.head(2)


# ## Read mutation data
# 
# This file contains mutation data (which mutations each sample contains). See the [online documentation](https://xenabrowser.net/datapages/?dataset=mc3.v0.2.8.PUBLIC.xena&host=https%3A%2F%2Fpancanatlas.xenahubs.net) for the `MC3` resource. For more information about mutation calling refer to [Ellrott et al. 2018](https://doi.org/10.1016/j.cels.2018.03.002 "Scalable Open Science Approach for Mutation Calling of Tumor Exomes Using Multiple Genomic Pipelines")
# 
# Note that duplicate mutation rows, which [occur](https://groups.google.com/d/msg/ucsc-cancer-genomics-browser/eg6nJOFSefw/Z0BM6pU9BAAJ "Message on the Xena Browser Google Group") for samples that were sequenced multiple times, are filtered.

# In[14]:


# Load chromosome-symbol to entrez_gene_id mapping
path = os.path.join('download', 'genes', 'chromosome-symbol-mapper.tsv')
chr_sym_map_df = pandas.read_table(path)
chr_sym_map_df.chromosome = 'chr' + chr_sym_map_df.chromosome
chr_sym_map_df.head(2)


# In[15]:


path = os.path.join('download', 'mc3.v0.2.8.PUBLIC.xena.tsv.gz')
snp_mutation_df = (
    pandas.read_table(path)
    .rename(columns={'sample': 'sample_id', 'gene': 'symbol'})
    .merge(chr_sym_map_df, on='symbol')
    .drop(['chromosome', 'symbol'], axis='columns')
    .drop_duplicates()
)
snp_mutation_df.head(2)


# In[16]:


# Number of samples with at least one mutation
samples_with_mutation_calls = sorted(set(snp_mutation_df.sample_id))
len(samples_with_mutation_calls)


# In[17]:


# Number of genes with at least one mutation
snp_mutation_df.entrez_gene_id.nunique()


# In[18]:


# Mutations counts by type
snp_mutation_df.effect.value_counts().reset_index()


# ### Convert SNP mutations to gene mutations
# 
# The next cell specifies which mutations to preserve as gene-affecting, which were chosen according to the red & blue [mutation effects in Xena](http://xena.ucsc.edu/how-we-characterize-mutations/).

# In[19]:


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


# In[20]:


# Mutations effects that were observed but nut included
set(snp_mutation_df.effect.unique()) - mutations


# In[21]:


gene_mutation_df = (snp_mutation_df
    .query("effect in @mutations")
    .groupby(['sample_id', 'entrez_gene_id'])
    .apply(len)
    .reset_index()
    .rename(columns={0: 'count'})
)

gene_mutation_df.head(2)


# In[22]:


# Create a sample (rows) by gene (columns) matrix of mutation status
gene_mutation_mat_df = (gene_mutation_df
    .pivot_table(index='sample_id',
                 columns='entrez_gene_id',
                 values='count',
                 fill_value=0)
    .reindex(samples_with_mutation_calls, fill_value=0)
    .astype(bool).astype(int)
)
gene_mutation_mat_df.columns = gene_mutation_mat_df.columns.astype(str)
gene_mutation_mat_df.shape


# In[23]:


'{:.2%} sample-gene pairs are mutated'.format(
    gene_mutation_mat_df.stack().mean())


# In[24]:


# Save complete mutation matrix
path = os.path.join('data', 'complete', 'mutation-matrix.tsv.bz2')
gene_mutation_mat_df.to_csv(path, sep='\t', compression='bz2')


# ## Read gene expression data
# 
# This file contains gene expression data. See the [online documentation](https://xenabrowser.net/datapages/?dataset=EB%2B%2BAdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.xena&host=https%3A%2F%2Fpancanatlas.xenahubs.net) for processing information.

# In[25]:


# Read the gene × sample dataset
path = os.path.join('download', 'EB++AdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.xena.tsv.gz')
expr_df = pandas.read_table(path, index_col=0)


# In[26]:


# Retrieve symbol to gene mapping for gene expression data
path = os.path.join('mapping', 'HiSeqV2-genes', 'HiSeqV2-gene-map.tsv')
gene_map_df = pandas.read_table(path, dtype='str')
symbol_to_entrez = dict(zip(gene_map_df.symbol, gene_map_df.entrez_gene_id))

# Check for unmapped symbols
unmapped_symbols = set(expr_df.index) - set(symbol_to_entrez)
unmapped_symbols


# In[27]:


# Process the dataset
expr_df = (expr_df
    # Drop genes with NA measurements
    # (See https://github.com/cognoma/cancer-data/issues/40)
    .dropna(axis='rows')
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

# Filter for valid Entrez Genes (This will filter the unmapped symbols)
expr_df = expr_df.loc[:, expr_df.columns.isin(gene_df.entrez_gene_id)]

expr_df.shape


# In[28]:


# Number of patients represented in the expression dataset
clinmat_df.query("sample_id in @expr_df.index").patient_id.nunique()


# In[29]:


# Peak at the data matrix
expr_df.iloc[:5, :5]


# In[30]:


# Save complete expression matrix
path = os.path.join('data', 'complete', 'expression-matrix.tsv.bz2')
expr_df.to_csv(path, sep='\t', float_format='%.3g', compression='bz2')


# ## Integrate expression and mutation data
# 
# Find samples with both mutation and expression data.
# 
# We assume that if a sample was not in the `MC3` data, it was not assayed for mutation ([more info](https://github.com/cognoma/cancer-data/issues/43#issuecomment-380957274)).

# In[31]:


sample_ids = sorted(clinmat_df.index & gene_mutation_mat_df.index & expr_df.index)
len(sample_ids)


# In[32]:


# Filter expression (x) and mutation (y) matrices for common samples
sample_df = clinmat_df.loc[sample_ids, :]
x_df = expr_df.loc[sample_ids, :]
y_df = gene_mutation_mat_df.loc[sample_ids, :]


# In[33]:


# Add a columnn for mutations per sample
sample_df['n_mutations'] = y_df.sum(axis='columns')


# In[34]:


x_gene_df = (
    gene_df.loc[x_df.columns, :]
    .assign(mean_expression=x_df.mean(axis='rows'))
    .assign(stdev_expression=x_df.std(axis='rows'))
)
path = os.path.join('data', 'expression-genes.tsv')
x_gene_df.to_csv(path, sep='\t', index=False, float_format='%.4g')
x_gene_df.head(2)


# In[35]:


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

# In[36]:


path = os.path.join('data', 'samples.tsv')
sample_df.to_csv(path, sep='\t', float_format='%.0f', index=False)


# In[37]:


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


# In[38]:


tsv_args = {'sep': '\t', 'float_format': '%.3g'}

for df, name in (x_df, 'expression-matrix'), (y_df, 'mutation-matrix'):

    # Save full dataset
    path = os.path.join('data', name + '.tsv.bz2')
    df.to_csv(path, **tsv_args, compression='bz2')
    
    # Save subsetted datasets
    for sample, nrows, ncols in ('small', 50, 15), ('all-samples', None, 15), ('all-genes', 50, None):
        path = os.path.join('data', 'subset', '{}-{}.tsv'.format(name, sample))
        subset_df(df, nrows=nrows, ncols=ncols).to_csv(path, **tsv_args)

