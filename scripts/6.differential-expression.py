
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

def get_diffex(subtype_df):
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
        ('entrez_gene_id', diffex_df.columns.astype(int)),
        ('patients', len(diffex_df)),
        ('tumor_mean', tumor_df.mean()),
        ('normal_mean', normal_df.mean()),
        ('mean_diff', diffex_df.mean()),
        ('t_stat', ttest.statistic),
        ('mlog10_p_value', -numpy.log10(ttest.pvalue)),
    ])
    return df

diffex_df = type_df.groupby('acronym').apply(get_diffex).reset_index('acronym')


# In[6]:

# Add gene symbols
path = os.path.join('data', 'genes.tsv')
gene_df = pandas.read_table(path)
gene_df = gene_df[['entrez_gene_id', 'symbol']]
diffex_df = gene_df.merge(diffex_df, how='right')


# In[7]:

diffex_df.head()


# In[8]:

path = os.path.join('data', 'complete', 'differential-expression.tsv.bz2')
diffex_df.to_csv(path, sep='\t', index=False, compression='bz2', float_format='%.4g')

