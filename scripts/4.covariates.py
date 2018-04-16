
# coding: utf-8

# # Creating a covariate dataset that encodes categorical sample variables

# The sample data possesses covariates as shown in <a href='https://github.com/cognoma/cancer-data/blob/master/3.explore-mutations.ipynb'>this notebook</a>. These may provide a spurious signal that a classifier accommodates for, and could confound attempts to pick out the actual signal that we desire. This notebook will create a file with encoded information on these covariates. Classifiers being implemented by the machine learning group can use this as additional data to train on.

# In[1]:


import os

import pandas as pd
import numpy as np


# Let's peak at the sample data.

# In[2]:


path = os.path.join('data', 'samples.tsv')

# Protect against int columns with missing values getting converted to floats
dtypes = {
    'dead': str,
    'recurred': str,
}
covariates_df = pd.read_table(path, index_col=0, dtype=dtypes)
covariates_df['recurred'] = covariates_df.days_recurrence_free.notnull().astype(int)
covariates_df.head(4).transpose()


# Let's get a sense of how much missing data there is.

# In[3]:


print('Total number of samples: {:,}'.format(len(covariates_df)))
print('Number of nulls in each column:')
covariates_df.isnull().sum(axis=0)


# In[4]:


# Specify which variables are categorical
categorical_variables = ['acronym', 'gender', 'dead', 'recurred']

# Number of categories per categorical variable
covariates_df[categorical_variables].apply(lambda x: x.nunique())


# In[5]:


continuous_variables = ['age_diagnosed', 'days_survived', 'days_recurrence_free', 'n_mutations']
covariates_df = covariates_df[categorical_variables + continuous_variables]
covariates_df.head(2)


# Before encoding, we're going to use the disease categories below. So let's store them for later use.

# Inspecting the head of the covariates DataFrame above, we see that two columns, namely <code>dead</code> and <code>recurred</code>, need some attention. We're going to encode categorical variables using panda's get_dummies. Since these columns are indicated by a 1 or 0, this will become the column header when encoded, as below.

# Let's rename the values in each of these so that they more accurately reflect the underlying variable.

# In[6]:


renamer = {
    'dead_0': 'alive',
    'dead_1': 'dead',
    'recurred_0': 'has_not_recurred',
    'recurred_1': 'has_recurred',
    'gender_Female': 'female',
    'gender_Male': 'male',
}
covariates_df = pd.get_dummies(covariates_df, columns=categorical_variables).rename(columns=renamer)
covariates_df.head(2)


# Now the column name more accurately reflects the underlying variable. The categorical values are encoded as numeric data that can be input to the types of classifiers that we have been using.
# 
# Another useful covariate will be the logarithm plus one function of the number mutations that was calculated in the aforementioned notebook.

# In[7]:


covariates_df['n_mutations_log1p'] = np.log1p(covariates_df.pop('n_mutations'))
covariates_df.head(2)


# Finally, let's save this to a <code>.tsv</code> file.

# In[8]:


path = os.path.join('data', 'covariates.tsv')
covariates_df.to_csv(path, sep='\t', float_format='%.5g')

