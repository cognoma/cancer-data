
# coding: utf-8

# # Download Cognoma genes

# In[1]:

import os
from urllib.request import urlretrieve


# In[2]:

repo = 'cognoma/genes'
commit = '721204091a96e55de6dcad165d6d8265e67e2a48'
base_url = 'https://github.com/{repo}/raw/{commit}/data/'.format(repo=repo, commit=commit)


# In[3]:

filenames = [
    'genes.tsv',
    'updater.tsv',
    'chromosome-symbol-mapper.tsv',
]

gene_dir = os.path.join('download', 'genes')

if not os.path.isdir(gene_dir):
    os.makedirs(gene_dir, exist_ok=True)
    
for name in filenames:
    url = base_url + name
    path = os.path.join(gene_dir, name)
    urlretrieve(url, path)

