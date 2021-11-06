#!/usr/bin/env python
# coding: utf-8

# In[2]:


import matplotlib.pyplot as plt
import pandas as pd
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from sklearn.model_selection import train_test_split


# In[3]:


Dataset = pd.read_csv('../../../normalizedCounts.csv')
Dataset_unnor = pd.read_csv('../../../UnnormalizedCounts.csv')


# In[4]:


list(Dataset)
# df_new = Dataset.rename(columns={'A': 'Col_1'}, index={'ONE': 'Row_1'})
df_new = Dataset.rename(columns={'Unnamed: 0': ''})


# In[5]:


# Save to file in csv format
df_new.to_csv('../../../NormalizedCountWithGeneColumnRenamed.csv')


# In[6]:


df_index = df_new.set_index('')
df_index
# Transpose Data
df_data = df_index.T
df_data


# In[7]:


df_data['labels'] = [
    'subacute ischemic stroke', 
    'subacute ischemic stroke', 
    'subacute ischemic stroke', 
    'subacute ischemic stroke', 
    'subacute ischemic stroke', 
    'acute ischemic stroke', 
    'acute ischemic stroke', 
    'acute ischemic stroke', 
    'acute ischemic stroke', 
    'acute ischemic stroke', 
    'healthy control', 
    'healthy control', 
    'healthy control', 
    'healthy control', 
    'healthy control'
]


# In[8]:


# Save new data

df_data.to_csv('../../../NormalizedCountWithGeneColumnLabels.csv')


# In[9]:


# df_data.describe()


# In[10]:


header = list(df_data.columns)
head_ = header.pop(60237)
len(header)


# In[11]:


df_data.loc[:,header].values


# In[12]:


# Separating out the features
df_x = df_data.loc[:, header].values
# Separating out the target
df_y = df_data.loc[:,['labels']].values
# Standardizing the features
x = StandardScaler().fit_transform(df_x)


# In[13]:


df_x_x = df_data.loc[:, df_data.isin([' ','NULL',0]).mean() < 1]


# In[14]:


df_x_x.shape


# In[15]:


df_x_x_zeros = df_data.loc[:, df_data.isin([' ','NULL',0]).mean() == 1]


# In[16]:


df_x_x_zeros


# In[17]:


df_data.isin([' ','NULL',0])['ENSG00000000005']
df_data.isin([' ','NULL',0])['ENSG00000000419']


# In[18]:


# Remove Label for analysis

df_without_label = df_x_x.drop(['labels'], axis=1)


# In[19]:


df_without_label.shape


# In[20]:


df3 = df_without_label[df_without_label.columns[df_without_label.mean() > 0.1]]


# In[21]:


df3.shape


# In[22]:


df3_less_than_0_1 = df_without_label[df_without_label.columns[df_without_label.mean() < 0.1]]


# In[23]:


df3_less_than_0_1


# In[24]:


# Get column data
df3['ENSG00000288640']

# get row data
df3.loc['GSM3483766']


# In[25]:


pca = PCA(n_components=3)
principalComponents = pca.fit_transform(df3)
principalComponents


# In[26]:


principal_df = pd.DataFrame(data = principalComponents
             , columns = ['p1', 'p2', 'p3'],
            index=df3.index.values)


# In[27]:


principal_df


# In[28]:


df3.index.values


# In[29]:


df_data[['labels']]


# In[30]:


final_df = pd.concat([principal_df, df_data[['labels']]], axis = 1)


# In[31]:


final_df


# In[32]:


fig = plt.figure(figsize = (8,8))
ax = fig.add_subplot(1,1,1) 
ax.set_xlabel('p1', fontsize = 15)
ax.set_ylabel('p2', fontsize = 15)
ax.set_title('2 component PCA', fontsize = 20)
targets = ['subacute ischemic stroke', 'acute ischemic stroke', 'healthy control']
colors = ['r', 'g', 'b']
for target, color in zip(targets,colors):
    indicesToKeep = final_df['labels'] == target
    ax.scatter(final_df.loc[indicesToKeep, 'p1']
               , final_df.loc[indicesToKeep, 'p2']
               , c = color
               , s = 50)
ax.legend(targets)
ax.grid()


# In[33]:


pca.explained_variance_ratio_


# In[34]:


print(pca.explained_variance_)
print(pca.explained_variance_ratio_)
print(pca.explained_variance_ratio_.cumsum())


# In[35]:


header_df3 = list(df3.columns)


# In[36]:


train_data, test_data, train_lbl, test_lbl = train_test_split( df3.loc[:,header_df3].values, 
                                                              final_df['labels'], 
                                                              test_size=1/7.0, 
                                                              random_state=42, shuffle=True)


# In[37]:


scaler = StandardScaler()
# Fit on training set only.
scaler.fit(train_data)
# Apply transform to both the training set and the test set.
train_data_ = scaler.transform(train_data)
test_data_ = scaler.transform(test_data)


# In[38]:


pca = PCA(.60)


# In[39]:


pca
pca.fit(train_data_)


# In[40]:


# Get predicted classes
pca.n_components_


# In[41]:


train_data_transform = pca.transform(train_data_)
test_data_transform = pca.transform(test_data_)


# In[42]:


from lazypredict.Supervised import LazyClassifier, LazyRegressor


# In[ ]:


reg = LazyRegressor(verbose=0,ignore_warnings=True,
custom_metric=None, predictions=True)
models,pred = reg.fit(train_data, test_data, train_lbl, test_lbl)


# In[ ]:





# In[ ]:





# In[ ]:
!pip install lazypredict




# %%
