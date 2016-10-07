
# coding: utf-8

# # Download Cognoma genes

# In[1]:

import os
from urllib.request import urlretrieve


# In[2]:

repo = 'dhimmel/genes'
commit = 'b64fcb4261005cf717d5d5ef4e03540e4a1f361e'
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

