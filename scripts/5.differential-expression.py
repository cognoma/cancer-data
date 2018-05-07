
# coding: utf-8

# ## Compute cancer-type-specific differential expression scores for each gene
# 
# Compares tumor to normal tissue in the same patient.

# In[1]:

import os

import pandas
import numpy
from scipy.stats import ttest_1samp


# In[2]:

path = os.path.join('data', 'complete', 'expression-matrix.tsv.bz2')
expr_df = pandas.read_table(path, index_col=0)


# In[3]:

path = os.path.join('data', 'complete', 'samples.tsv')
sample_df = (
    pandas.read_table(path)
    # Filter for samples with expression
    .query("sample_id in @expr_df.index")
)
patient_df = sample_df[['patient_id', 'acronym']].drop_duplicates()
sample_df.head(2)


# In[4]:

type_df = sample_df.pivot('patient_id', 'sample_type', values='sample_id')
type_df = type_df[['Primary Tumor', 'Solid Tissue Normal']]
# Filter for paired samples
type_df = type_df[type_df.isnull().sum(axis='columns') == 0]
type_df = type_df.reset_index().merge(patient_df)
type_df.head(2)


# In[5]:

def get_diffex(subtype_df, expr_df):
    """
    For each gene, compute differential expression between paired tumor and normal tissue.
    """
    tumor_df = expr_df.loc[list(subtype_df['Primary Tumor']), :]
    normal_df = expr_df.loc[list(subtype_df['Solid Tissue Normal']), :]
    for df in tumor_df, normal_df:
        df.index = subtype_df.index
    
    diffex_df = tumor_df - normal_df
    ttest = ttest_1samp(diffex_df, popmean=0, axis=0)

    df = pandas.DataFrame.from_items([
        ('entrez_gene_id', diffex_df.columns),
        ('patients', len(diffex_df)),
        ('tumor_mean', tumor_df.mean()),
        ('normal_mean', normal_df.mean()),
        ('mean_diff', diffex_df.mean()),
        ('t_stat', ttest.statistic),
        ('mlog10_p_value', -numpy.log10(ttest.pvalue)),
    ])
    return df


# In[6]:

diffex_df = (type_df
    .groupby('acronym')
    .apply(get_diffex, expr_df=expr_df)
    .reset_index('acronym')
    .query("patients >= 5")
)

diffex_df.entrez_gene_id = diffex_df.entrez_gene_id.astype(int)


# In[7]:

# Add gene symbols
path = os.path.join('data', 'expression-genes.tsv')
gene_df = pandas.read_table(path, low_memory=False)
gene_df = gene_df[['entrez_gene_id', 'symbol']]
len(gene_df)


# In[8]:

diffex_df = diffex_df.merge(gene_df, how='left')
diffex_df.tail()


# In[9]:

path = os.path.join('data', 'complete', 'differential-expression.tsv.bz2')
diffex_df.to_csv(path, sep='\t', index=False, compression='bz2', float_format='%.4g')


# In[10]:

# Patients with paired samples per disease
path = os.path.join('download', 'diseases.tsv')
acronym_df = pandas.read_table(path)

(type_df.acronym
    .value_counts().rename('patients')
    .reset_index().rename(columns={'index': 'acronym'})
    .merge(acronym_df)
    .sort_values('patients', ascending=False)
)


# # Reduce expression dimensionality with NMF

# In[11]:

from sklearn.decomposition import NMF
import matplotlib.pyplot as plt
import seaborn

get_ipython().magic('matplotlib inline')


# In[12]:

nmf = NMF(n_components=100, random_state=0)
expr_nmf_df = nmf.fit_transform(expr_df)
expr_nmf_df = pandas.DataFrame(expr_nmf_df, index=expr_df.index)
expr_nmf_df.iloc[:5, :5]


# In[13]:

diffex_nmf_df = (type_df
    .groupby('acronym')
    .apply(get_diffex, expr_df=expr_nmf_df)
    .reset_index('acronym')
    .rename(columns={'entrez_gene_id': 'component'})
    .query("patients >= 5")
)
diffex_nmf_df.head(2)


# In[14]:

# Acronyms at https://github.com/cognoma/cancer-data/blob/master/download/diseases.tsv
plot_df = diffex_nmf_df.pivot(index='acronym', columns='component', values='t_stat').fillna(0)
grid = seaborn.clustermap(plot_df, metric='correlation', figsize=(9, 6))
grid.ax_heatmap.tick_params(labelbottom='off')
_ = plt.setp(grid.ax_heatmap.yaxis.get_majorticklabels(), rotation=0)

