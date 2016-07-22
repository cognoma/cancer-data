
# coding: utf-8

# # Download TCGA Pan-Cancer Datasets from the UCSC Xena Browser
# 
# This notebook downloads TCGA datasets for Project Cognoma. The file contents (text) remains unmodified, but files are given extensions and bzip2 compressed.
# 
# [See here](https://genome-cancer.soe.ucsc.edu/proj/site/xena/datapages/?cohort=TCGA%20Pan-Cancer%20%28PANCAN%29 "Xena: cohort: TCGA Pan-Cancer (PANCAN)") for all TCGA Pan-Cancer datasets on Xena.

# In[1]:

import os
import bz2
from urllib.request import urlretrieve


# In[2]:

def bzip2_compress(path, keep=False):
    """
    Compress a file using bzip2 compression.
    Designed to mirror the functionality of the
    `bzip2 --compress $PATH` shell command.
    `keep` specifies whether to remove the uncompressed file.
    """
    with open(path, 'rb') as reader, bz2.open(path + '.bz2', 'wb') as writer:
        writer.writelines(reader)
    if not keep:
        os.remove(path)


# Documentation for the TCGA Pan-Cancer files from the Xena browser:
# 
# + [`HiSeqV2`](https://genome-cancer.soe.ucsc.edu/proj/site/xena/datapages/?dataset=TCGA.PANCAN.sampleMap/HiSeqV2&host=https://tcga.xenahubs.net)
# + [`PANCAN_clinicalMatrix`](https://genome-cancer.soe.ucsc.edu/proj/site/xena/datapages/?dataset=TCGA.PANCAN.sampleMap/PANCAN_clinicalMatrix&host=https://tcga.xenahubs.net)
# + [`PANCAN_mutation`](https://genome-cancer.soe.ucsc.edu/proj/site/xena/datapages/?dataset=TCGA.PANCAN.sampleMap/PANCAN_mutation&host=https://tcga.xenahubs.net)

# In[3]:

base_url = 'https://tcga.xenahubs.net/download/TCGA.PANCAN.sampleMap/'

names = [
    'PANCAN_clinicalMatrix',
    'PANCAN_mutation',
    'HiSeqV2',
]


# In[4]:

# Download metadata
for name in names:
    url = base_url + name + '.json'
    path = os.path.join('download', name + '.json')
    urlretrieve(url, path)


# In[5]:

# Download data
for name in names:
    url = base_url + name
    path = os.path.join('download', name + '.tsv')
    urlretrieve(url, path)
    bzip2_compress(path)

