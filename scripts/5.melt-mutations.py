
# coding: utf-8

# # Melt mutation data

# This notebook melts the mutation matrix into a tabular format, for use by the <code>core-service</code> repository.
# 
# Imports...

# In[1]:

import pandas as pd
import os


# Load mutation matrix...

# In[2]:

get_ipython().run_cell_magic('time', '', "mutation_path = os.path.join('data', 'mutation-matrix.tsv.bz2')\nmutation_df = pd.read_table(mutation_path)")


# Grab a small subset of mutated genes to explain the approach...

# In[3]:

mutation_subset_df = mutation_df[mutation_df['1']==1].head().loc[:, :'12']
mutation_subset_df


# <a href='http://pandas.pydata.org/pandas-docs/stable/generated/pandas.melt.html'>Melt</a> the mutation data frame. This will convert to a table of <code>sample</code>, <code>entrez_gene_id</code>, <code>mutation_status</code> pairs. Then, <code>mutation_status</code> is filtered to only contain mutated sample/gene pairs. Note that it correctly picks out the seven mutated sample/gene pairs.

# In[4]:

melted_mutation_subset_df = pd.melt(mutation_subset_df, id_vars='sample_id', var_name='entrez_gene_id')
melted_mutation_subset_df[melted_mutation_subset_df.value==1].drop('value', axis=1)


# Finally, apply this to the full mutation data, and write to a <code>.tsv</code> file.

# In[5]:

melted_mutation_path = os.path.join('data', 'melted-mutations.tsv')
melted_mutation_df = pd.melt(mutation_df, id_vars='sample_id', var_name='entrez_gene_id')
mutations = melted_mutation_df.value==1
melted_mutation_df[mutations].drop('value', axis=1).to_csv(melted_mutation_path, sep='\t', index=False)

