
# coding: utf-8

# # Creating a covariate dataset that encodes categorical sample variables

# The sample data possesses covariates as shown in <a href='https://github.com/cognoma/cancer-data/blob/master/3.explore-mutations.ipynb'>this notebook</a>. These may provide a spurious signal that a classifier accommodates for, and could confound attempts to pick out the actual signal that we desire. This notebook will create a file with encoded information on these covariates. Classifiers being implemented by the machine learning group can use this as additional data to train on.

# In[1]:

import os
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
get_ipython().magic('matplotlib inline')


# Let's peak at the sample data.

# In[2]:

path = os.path.join('data', 'samples.tsv')
covariates_df = pd.read_table(path, index_col=0)
covariates_df.head(10)


# Let's get a sense of how much missing data there is.

# In[3]:

print('Total number of samples: %d' % len(covariates_df))
print('Number of nulls in each column:')
covariates_df.isnull().sum(axis=0)


# Now count the number of types of each column. We'll use Python's set, and therefore must exclude null values.

# In[4]:

categorical_variables = ['patient_id', 'sample_type', 'disease', 'organ_of_origin', 'gender', 'dead', 'recurred']
for variable in categorical_variables:
    not_null = covariates_df[variable].notnull()
    count = len(set(covariates_df.loc[not_null, variable]))
    print('Number of types in category %s: %d' % (variable, count))


# There are no missing values from <code>sample_type</code> and only one possible value, so it is redundant and we'll exclude it. Also, <code>patient_id</code> has as many classes as there are samples, so this is unique and will be excluded as well.

# In[5]:

covariates_df.drop(['patient_id', 'sample_type'], axis=1, inplace=True)
del categorical_variables[:2]
print('Remaining variables to encode: %s' % ', '.join(categorical_variables))


# Before encoding, we're going to use the disease categories below. So let's store them for later use.

# In[6]:

diseases = covariates_df.disease


# Inspecting the head of the covariates DataFrame above, we see that two columns, namely <code>dead</code> and <code>recurred</code>, need some attention. We're going to encode categorical variables using panda's get_dummies. Since these columns are indicated by a 1 or 0, this will become the column header when encoded, as below.

# In[7]:

pd.get_dummies(covariates_df, columns=categorical_variables).head(2)


# Let's rename the values in each of these so that they more accurately reflect the underlying variable.

# In[8]:

renamer = {'dead_0.0':'alive',
           'dead_1.0':'dead',
           'recurred_0.0':'has_not_recurred',
           'recurred_1.0': 'has_recurred'}
covariates_df = pd.get_dummies(covariates_df, columns=categorical_variables).rename(columns=renamer)
covariates_df.head(2)


# Now the column name more accurately reflects the underlying variable. The categorical values are encoded as numeric data that can be input to the types of classifiers that we have been using.
# 
# Another useful covariate will be the logarithm plus one function of the number mutations that was calculated in the aforementioned notebook.

# In[9]:

path = os.path.join('data', 'mutation-matrix.tsv.bz2')
mutation_df = pd.read_table(path, index_col=0)


# In[10]:

mutation_df.head(2)


# In[11]:

number_of_mutations = mutation_df.sum(axis='columns')
number_of_mutations.head(2)


# Let's add this as an additional column to our covariates.

# In[12]:

covariates_df['n_mutations_log1p'] = np.log1p(number_of_mutations)
covariates_df.head()


# Finally, let's save this to a <code>.tsv</code> file.

# In[13]:

path = os.path.join('data', 'covariates.tsv')
covariates_df.to_csv(path, sep='\t', float_format='%.5g')

