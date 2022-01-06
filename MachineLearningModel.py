import sys
import matplotlib.pyplot as plt
import pandas as pd
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from sklearn.model_selection import train_test_split
import seaborn as sns
import os
from matplotlib.backends.backend_pdf import PdfPages
def machineLearning(
    file="../../../normalizedCounts.csv", 
    path_to_save="machineLearning", 
    component_n=3, 
    condition="SUB,NC", 
    sample_count="5, 5",
    graphFileName="graph"
):
    pp = PdfPages('{}/{}.pdf'.format(path_to_save, graphFileName))
    component_n = int(component_n)
    # Check whether the specified path exists or not
    isExist = os.path.exists(path_to_save)

    if not isExist:
        # Create a new directory because it does not exist 
        os.makedirs(path_to_save)
    print("{} is created!".format(path_to_save))
    condition = condition.split(',')
    sample_count = sample_count.split(',')
    if len(sample_count) != len(condition):
        print("parameter slice_group, condition and sample_count must have the same length. ")
    Dataset = pd.read_csv(file)
    list(Dataset)
    df_new = Dataset.rename(columns={'Unnamed: 0': ''})

    # Save to file in csv format
    df_new.to_csv('{}/NormalizedCountWithGeneColumnRenamed.csv'.format(path_to_save))

    df_index = df_new.set_index('')
    df_index
    # Transpose Data
    df_data = df_index.T
    df_data

    dictionary = dict(zip(condition, sample_count))
    category_dup = []
    for key, val in dictionary.items():
        category_dup += [key]*int(val)
    df_data['labels'] = category_dup
    print(df_data['labels'])

    # Save new data

    df_data.to_csv('./{}/NormalizedCountWithGeneColumnLabels.csv'.format(path_to_save))



    # df_data.describe()


    header = list(df_data.columns)
    head_ = header.pop(60237)
    len(header)


    df_data.loc[:,header].values

    # Separating out the features
    df_x = df_data.loc[:, header].values
    # Separating out the target
    df_y = df_data.loc[:,['labels']].values
    # Standardizing the features
    x = StandardScaler().fit_transform(df_x)



    df_x_x = df_data.loc[:, df_data.isin([' ','NULL',0]).mean() < 1]


    df_x_x.shape

    df_x_x_zeros = df_data.loc[:, df_data.isin([' ','NULL',0]).mean() == 1]

    df_x_x_zeros

    # df_data.isin([' ','NULL',0])['ENSG00000000005']
    # df_data.isin([' ','NULL',0])['ENSG00000000419']


    # Remove Label for analysis

    df_without_label = df_x_x.drop(['labels'], axis=1)

    df_without_label.shape


    df3 = df_without_label[df_without_label.columns[df_without_label.mean() > 0.1]]

    df3.shape

    df3_less_than_0_1 = df_without_label[df_without_label.columns[df_without_label.mean() < 0.1]]

    df3_less_than_0_1



    # Get column data
    # df3['ENSG00000288640']

    # # get row data
    # df3.loc['GSM3483766']

    pca = PCA(n_components=component_n)
    principalComponents = pca.fit_transform(df3)
    principalComponents

    principal_df = pd.DataFrame(data = principalComponents
                , columns = ['P'+ str(i) for i in range(1, component_n+1)],
                index=df3.index.values)


    principal_df


    df3.index.values


    df_data[['labels']]


    final_df = pd.concat([principal_df, df_data[['labels']]], axis = 1)

    final_df

    fig = plt.figure(figsize = (8,8))
    ax = fig.add_subplot(1,1,1) 
    ax.set_xlabel('p1', fontsize = 15)
    ax.set_ylabel('p2', fontsize = 15)
    ax.set_title('2 component PCA', fontsize = 20)
    targets = condition
    colors = ['r', 'g', 'b']
    for target, color in zip(targets,colors):
        indicesToKeep = final_df['labels'] == target
        ax.scatter(final_df.loc[indicesToKeep, 'P1']
                , final_df.loc[indicesToKeep, 'P2']
                , c = color
                , s = 50)
    ax.legend(targets)
    ax.grid()

    pp.savefig(fig)


    pca.explained_variance_ratio_


    print(pca.explained_variance_)
    print(pca.explained_variance_ratio_)
    print(pca.explained_variance_ratio_.cumsum())

    header_df3 = list(df3.columns)

    train_data, test_data, train_lbl, test_lbl = train_test_split( df3.loc[:,header_df3].values, 
                                                                final_df['labels'], 
                                                                test_size=1/7.0, 
                                                                random_state=42, shuffle=True)


    scaler = StandardScaler()
    # Fit on training set only.
    scaler.fit(train_data)
    # Apply transform to both the training set and the test set.
    train_data_ = scaler.transform(train_data)
    test_data_ = scaler.transform(test_data)

    pca = PCA(.60)

    pca
    pca.fit(train_data_)


    # Get predicted classes
    pca.n_components_

    train_data_transform = pca.transform(train_data_)
    test_data_transform = pca.transform(test_data_)


    from lazypredict.Supervised import LazyClassifier, LazyRegressor


    # reg = LazyClassifier(predictions=True)
    # models,pred = reg.fit(train_data, test_data, train_lbl, test_lbl)

    clf = LazyClassifier(verbose=0,ignore_warnings=True, custom_metric=None, predictions=True)
    regr=LazyRegressor(verbose=0,predictions=True)

    print("Here we are")
    print(train_lbl.shape)
    print(test_lbl.shape)
    models_train,predictions_train = clf.fit(train_data_transform, train_data_transform, train_lbl, train_lbl)
    models_test,predictions_test = regr.fit(train_data_transform, test_data_transform, train_lbl, test_lbl)


    print(predictions_train)
    print(models_test)
    print(predictions_test)

    pd.DataFrame(train_lbl).to_csv('./{}/label_train.csv'.format(path_to_save))
    pd.DataFrame(test_lbl).to_csv('./{}/label_test.csv'.format(path_to_save))
    predictions_train.to_csv('./{}/predictions_train.csv'.format(path_to_save))
    models_test.to_csv('./{}/models_test.csv'.format(path_to_save))
    predictions_test.to_csv('./{}/predictions_test.csv'.format(path_to_save))
    models_train.to_csv('./{}/modelTained.csv'.format(path_to_save))
    
    fig1 = plt.figure(figsize=(30, 20))
    sns.set_theme(style="whitegrid")
    ax = sns.barplot(x=models_train.index, y="Accuracy", data=models_train)
    plt.xticks(rotation=90)
    plt.savefig( '{}/{}.png'.format(path_to_save, 'AccuracyPlotH.png'), bbox_inches='tight')
    pp.savefig(fig1)


    fig2 = plt.figure(figsize=(20, 30))
    sns.set_theme(style="whitegrid")
    ax = sns.barplot(y=models_train.index, x="Accuracy", data=models_train)
    plt.savefig('{}/{}.png'.format(path_to_save, 'AccuracyPlot.png'), bbox_inches='tight')
    pp.savefig(fig2)


    fig3 = plt.figure(figsize=(20, 30))
    sns.set_theme(style="whitegrid")
    ax = sns.barplot(x=models_train.index, y="Time Taken", data=models_train)
    plt.xticks(rotation=90)
    plt.savefig('{}/{}.png'.format(path_to_save, 'TimePlot.png'), bbox_inches='tight')

    pp.savefig(fig3)



    # plt.figure(figsize=(10, 5))
    # sns.set_theme(style="whitegrid")
    # ax = sns.barplot(x=models_train.index, y="R-Squared", data=models_train)
    # ax.set(ylim=(0, 1))
    # plt.xticks(rotation=90)


    # plt.figure(figsize=(5, 10))
    # sns.set_theme(style="whitegrid")
    # ax = sns.barplot(y=models_train.index, x="R-Squared", data=models_train)
    # ax.set(xlim=(0, 1))

    # plt.savefig("{}/{}{}".format(path_to_save, graphFileName,'.pdf'), 
    #     dpi=None, 
    #     facecolor='w', 
    #     edgecolor='w',
    #     orientation='portrait', 
    #     papertype=None, 
    #     format=None,
    #     transparent=False, 
    #     bbox_inches=None, 
    #     pad_inches=0.1,
    #     frameon=None, 
    #     metadata=None)
    pp.close()



# machineLearning(
#     file="../../../UnnormalizedCounts.csv", 
#     path_to_save="machineLearning", 
#     component_n=3, 
#     condition="AIS, SUB,NC", 
#     sample_count="5, 5, 5",
#     graphFileName="graphs"
# )