
# coding: utf-8

# # Export cancer-data to JSON for frontend

# In[1]:

import os
import json
import gzip

import pandas


# ## gene_to_mutated_samples

# In[2]:

path = os.path.join('data', 'mutation-matrix.tsv.bz2')
mutation_mat_df = pandas.read_table(path, index_col=0)


# In[3]:

gene_to_mutated_samples = dict()
for entrez_gene_id, series in mutation_mat_df.iteritems():
    gene_to_mutated_samples[int(entrez_gene_id)] = list(series.index[series == 1])


# In[4]:

path = os.path.join('data', 'json', 'gene_to_mutated_samples.json')
with open(path, 'w') as write_file:
    json.dump(gene_to_mutated_samples, write_file, indent=2, sort_keys=True)
'{:.2f} MB'.format(1e-6 * os.path.getsize(path))


# In[5]:

path = os.path.join('data', 'json', 'gene_to_mutated_samples.json.gz')
with gzip.open(path, 'wt') as write_file:
    json.dump(gene_to_mutated_samples, write_file, indent=2, sort_keys=True)
'{:.2f} MB'.format(1e-6 * os.path.getsize(path))


# ## disease_to_samples

# In[6]:

path = os.path.join('data', 'samples.tsv')
sample_df = pandas.read_table(path)
sample_df.head(2)


# In[7]:

disease_to_samples = {acronym: sorted(sample_ids) for acronym, sample_ids in sample_df.groupby('acronym').sample_id}


# In[8]:

path = os.path.join('data', 'json', 'disease_to_samples.json')
with open(path, 'w') as write_file:
    json.dump(disease_to_samples, write_file, indent=2, sort_keys=True)
'{:.2f} MB'.format(1e-6 * os.path.getsize(path))

