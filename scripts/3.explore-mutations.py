
# coding: utf-8

# # Exploratory data analysis of TCGA mutation data

# In[1]:

import os

import numpy
import pandas
import seaborn

get_ipython().magic('matplotlib inline')


# ## Read TCGA datasets

# In[2]:

path = os.path.join('data', 'mutation-matrix.tsv.bz2')
mutation_df = pandas.read_table(path, index_col=0)
mutation_df.columns.name = 'entrez_gene_id'
mutation_df.shape


# In[3]:

path = os.path.join('data', 'samples.tsv')
sample_df = pandas.read_table(path)
sample_df.head(2)


# ## Distribution of mutations counts for genes

# In[4]:

gene_df = mutation_df.sum(axis='rows').rename('n_mutations').reset_index()
gene_df['n_mutations_log1p'] = numpy.log1p(gene_df.n_mutations)
gene_df.head(2)


# In[5]:

ax = seaborn.distplot(gene_df.n_mutations_log1p)
xticks = ax.get_xticks()
xticklabels = numpy.expm1(xticks).round().astype(int)
axis_texts = ax.set_xticklabels(xticklabels)


# In[6]:

sum(gene_df.n_mutations == 0)


# ## Distribution of mutations counts for samples

# In[7]:

sample_df = sample_df.merge(
    mutation_df.sum(axis='columns').rename('n_mutations').reset_index()
)
sample_df['n_mutations_log1p'] = numpy.log1p(sample_df.n_mutations)
sample_df.head(2)


# In[8]:

# Mutations per sample
ax = seaborn.distplot(sample_df.n_mutations_log1p)
xticks = ax.get_xticks()
xticklabels = numpy.expm1(xticks).round().astype(int)
axis_texts = ax.set_xticklabels(xticklabels)


# ## Diagnosis age versus mutation count for samples

# In[9]:

grid = seaborn.jointplot('n_mutations_log1p', 'age_diagnosed', data=sample_df, kind='hex')
xticks = grid.ax_marg_x.get_xticks()
xticklabels = numpy.expm1(xticks).round().astype(int)
axis_texts = grid.ax_marg_x.set_xticklabels(xticklabels)


# ## Mutation frequency by disease

# In[10]:

genes = mutation_df.columns.tolist()
verbose_mutation_df = sample_df.merge(mutation_df.reset_index())
mutation_freq_df = verbose_mutation_df.groupby('disease').apply(lambda df: df[genes].mean(axis='rows')).assign(
    n_mutations = verbose_mutation_df.groupby('disease').apply(len)
)
mutation_freq_df.iloc[:3, :3]


# In[11]:

gene_subset = {
    '7157': 'TP53', # tumor protein p53
    '7428': 'VHL', # von Hippel-Lindau tumor suppressor
    '29126': 'CD274', # CD274 molecule
}
plot_df = (mutation_freq_df
    .query("n_mutations > 100")
    [list(gene_subset)]
    .rename(columns=gene_subset)
)
ax = seaborn.heatmap(plot_df)

