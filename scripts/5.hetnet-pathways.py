
# coding: utf-8

# # Create a sample-by-pathway dataset

# This notebook attempts to build a dataset using the existing <code>mutation-matrix</code> dataset and the <a href='http://neo4j.het.io/browser/'>hetnet</a> that was introduced at the first meetup.
# 
# Imports...

# In[1]:

from neo4j.v1 import GraphDatabase
import pandas as pd
import numpy as np
import os
from ipywidgets import FloatProgress
from IPython.display import display


# Query the hetnet...

# In[2]:

query = '''
MATCH path = (gene:Gene)
  -[:PARTICIPATES_GpPW]-(pathway:Pathway)
RETURN gene.name as gene_name,
  gene.identifier as gene_id, 
  pathway.name as pathway_name,
  pathway.identifier as pathway_id
'''

driver = GraphDatabase.driver("bolt://neo4j.het.io")
hetnet_results = pd.DataFrame()
with driver.session() as session:
    result = session.run(query)
    hetnet_results = ( pd.DataFrame((x.values() for x in result), columns=result.keys())
                       .sort_values(by='pathway_id') )


# In[3]:

print('There were %d gene-pathway interactions in the hetnet query.' % len(hetnet_results))


# In[4]:

hetnet_results.head()


# In[5]:

path = os.path.join('data', 'mutation-matrix.tsv.bz2')
mutation_df = pd.read_table(path, index_col=0)
mutation_df.head()


# Note that columns of <code>mutation_df</code> are strings. After a quick search I couldn't find out how to fix this, so be warned that code below works around this.

# In[6]:

type(mutation_df.columns[0])


# Display some information about the two datasets...

# In[7]:

cognoma_genes = set([int(gene_id) for gene_id in mutation_df.columns])
print('There were %d genes found in the Cognoma dataset.' % len(cognoma_genes))

hetnet_genes = set(hetnet_results['gene_id'])
print('There were %d genes found in the hetnet pathway query.' % len(hetnet_genes))

genes_in_both = cognoma_genes.intersection(hetnet_genes)
print('There were %d genes found in both the Cognoma dataset and the hetnet query.' % len(genes_in_both))

pathways = set(hetnet_results['pathway_id'])
print('There were %d pathways found in the hetnet pathway query.' % len(pathways))

number_missing = len(cognoma_genes-hetnet_genes)
print('There were %d genes in the Cognoma dataset that were not found in the hetnet pathway query.' % number_missing)

number_missing = len(hetnet_genes-cognoma_genes)
print('There were %d genes in the hetnet pathway query that were not found in the Cognoma dataset.' % number_missing)


# Initialize a samples by pathways data frame...

# In[8]:

number_of_samples = len(mutation_df)
number_of_pathways = len(pathways) 
sample_pathway_df = pd.DataFrame(np.zeros((number_of_samples, number_of_pathways), dtype=np.int),
                                 index=mutation_df.index,
                                 columns=pathways)


# Now populate this data frame. This is a slow Python loop, hence the progress bar. It takes a few minutes on my laptop. The idea is to loop over all gene-pathway interactions in the hetnet query. If the gene is in the Cognoma dataset, we grab the pathway id in that gene-pathway interaction. We look at Cognoma samples where that gene is labeled 1, i.e., at Cognoma samples that have a mutation in that gene, and grab the corresponding indices. Then, in the pathway matrix all samples get the associated pathway tagged as a 1, since they have a mutated gene that participates in that pathway.

# In[9]:

i = 0
progress_bar = FloatProgress(min=0, max=len(hetnet_results))
display(progress_bar)
for _, row in hetnet_results.iterrows():
    gene_id = row['gene_id']
    if gene_id in genes_in_both:
        pathway_id = row['pathway_id']
        affected_samples = mutation_df.loc[:, str(gene_id)] == 1
        sample_pathway_df.loc[affected_samples, pathway_id] = 1
    i += 1
    progress_bar.value = i
sample_pathway_df.head()


# Finally, we write to disk. The raw file is about 26MB, so we use bz2 compression. The file is no longer tracked due to <code>data/.gitignore</code>.

# In[10]:

path = os.path.join('data','pathways.tsv.bz2')
sample_pathway_df.to_csv(path, sep='\t', compression='bz2')

