
# coding: utf-8

# # Export gene information

# In[1]:

import os
import collections

import pandas


# ## Compute local gene information

# In[2]:

path = os.path.join('data', 'mutation-matrix.tsv.bz2')
mutation_df = pandas.read_table(path, index_col=0)

path = os.path.join('data', 'expression-matrix.tsv.bz2')
expr_df = pandas.read_table(path, index_col=0)


# In[3]:

mutation_summary_df = pandas.DataFrame.from_items([
    ('entrez_gene_id', mutation_df.columns.astype(int)),
    ('n_mutations', mutation_df.sum(axis='rows')),
    ('mutation_frequency', mutation_df.mean(axis='rows')),
])

expression_summary_df = pandas.DataFrame.from_items([
    ('entrez_gene_id', expr_df.columns.astype(int)),
    ('mean_expression', expr_df.mean(axis='rows')),
])

summary_df = mutation_summary_df.merge(expression_summary_df, how='outer')
summary_df['mutation'] = summary_df.n_mutations.notnull().astype(int)
summary_df['expression'] = summary_df.mean_expression.notnull().astype(int)
summary_df.head(2)


# ## Retrieve Entrez Gene

# In[4]:

renamer = collections.OrderedDict([
    ('GeneID', 'entrez_gene_id'),
    ('Symbol', 'symbol'),
    ('description', 'description'),
    ('chromosome', 'chromosome'),
    ('type_of_gene', 'gene_type'),
    ('Synonyms', 'synonyms'),
    ('Other_designations', 'aliases'),
])

url = 'ftp://ftp.ncbi.nih.gov/gene/DATA/GENE_INFO/Mammalia/Homo_sapiens.gene_info.gz'
entrez_df = (pandas.read_table(url, compression='gzip')
    .rename(columns=renamer)
    [list(renamer.values())]
)
entrez_df.head()


# ## Combine local information with Entrez Gene

# In[5]:

combined_df = entrez_df.merge(summary_df, how='right').sort_values('entrez_gene_id')
len(combined_df)


# In[6]:

# Missing value per column
combined_df.isnull().sum()


# In[7]:

combined_df.head()


# In[8]:

# Write to TSV
path = os.path.join('data', 'genes.tsv')
combined_df.to_csv(path, sep='\t', index=False, float_format='%.4g')

