
# coding: utf-8

# # Download TCGA PanCanAtlas Datasets from the UCSC Xena Browser
# 
# This notebook downloads TCGA datasets for Project Cognoma. The file contents (text) remains unmodified, but files are given extensions and bzip2 compressed.
# 
# [See here](https://xenabrowser.net/datapages/?hub=https://pancanatlas.xenahubs.net:443 "Xena: cohort: TCGA PanCanAtlas") for all TCGA PanCanAtlas datasets on Xena.

# In[1]:


import os
from urllib.request import urlretrieve


# Documentation for the TCGA Pan-Cancer files from the Xena browser:
# 
# + Gene Expression: [`EB++AdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.xena`](https://xenabrowser.net/datapages/?dataset=EB%2B%2BAdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.xena&host=https%3A%2F%2Fpancanatlas.xenahubs.net)
# + Clinical Covariates: [`Survival_SupplementalTable_S1_20171025_xena_sp`](https://xenabrowser.net/datapages/?dataset=Survival_SupplementalTable_S1_20171025_xena_sp&host=https%3A%2F%2Fpancanatlas.xenahubs.net)
# + Mutation Data: [`mc3.v0.2.8.PUBLIC.xena.gz`](https://xenabrowser.net/datapages/?dataset=mc3.v0.2.8.PUBLIC.xena&host=https%3A%2F%2Fpancanatlas.xenahubs.net)

# In[2]:


base_url = 'https://pancanatlas.xenahubs.net/download/'

names = [
    'Survival_SupplementalTable_S1_20171025_xena_sp',
    'mc3.v0.2.8.PUBLIC.xena',
    'EB++AdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.xena',
]


# In[3]:


# Download metadata
for name in names:
    url = base_url + name + '.json'
    path = os.path.join('download', name + '.json')
    urlretrieve(url, path)


# In[4]:


# Download data
for name in names:
    url = base_url + name + '.gz'
    path = os.path.join('download', name + '.tsv.gz')
    urlretrieve(url, path)

