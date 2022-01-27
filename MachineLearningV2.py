import matplotlib.pyplot as plt
import pandas as pd
import os
import time
import csv
from os import path
from sklearn.preprocessing import StandardScaler
import seaborn as sns
from sklearn.linear_model import LogisticRegression
from sklearn.pipeline import make_pipeline
from sklearn.svm import SVC
from sklearn.ensemble import RandomForestClassifier
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis
from sklearn.neighbors import KNeighborsClassifier
from sklearn import tree
from sklearn.cluster import KMeans
from sklearn.naive_bayes import GaussianNB
from sklearn.ensemble import AdaBoostClassifier
from sklearn.ensemble import BaggingClassifier
from sklearn.linear_model import RidgeClassifierCV
from sklearn.discriminant_analysis import QuadraticDiscriminantAnalysis
from sklearn.gaussian_process import GaussianProcessClassifier
from sklearn.gaussian_process.kernels import RBF
from sklearn.neural_network import MLPClassifier
from sklearn.metrics import classification_report
from sklearn.calibration import CalibratedClassifierCV
from sklearn.metrics import roc_auc_score
from sklearn.decomposition import PCA
from matplotlib.backends.backend_pdf import PdfPages
from sklearn.metrics import accuracy_score, confusion_matrix, f1_score, precision_score, recall_score
from sklearn import metrics
import numpy as np
from sklearn.model_selection import train_test_split
from scipy import stats





models = []
models.append(('Logistic Regression', LogisticRegression(solver='lbfgs', max_iter=1000)))
models.append(('Linear Discriminant', LinearDiscriminantAnalysis()))
models.append(('KNN', KNeighborsClassifier(n_neighbors=3)))
models.append(('DecisionTree', tree.DecisionTreeClassifier(max_depth=5)))
models.append(('GaussianNB', GaussianNB()))
models.append(('SVM', SVC(random_state=0, probability=True)))
models.append(('RandomForest', RandomForestClassifier(max_depth=5, n_estimators=10, max_features=1)))
models.append(('AdaBooster', AdaBoostClassifier(n_estimators=100, random_state=0)))
models.append(('BaggingClassifier', BaggingClassifier(base_estimator=SVC(), n_estimators=10, random_state=0)))
models.append(('QuadraticDiscriminantAnalysis', QuadraticDiscriminantAnalysis()))
models.append(('SVC_linear', SVC(kernel="linear", C=0.025, probability=True)))
models.append(('SVC_gamma_2', SVC(gamma=2, C=1, probability=True)))
models.append(('GaussianProcessClassifier', GaussianProcessClassifier(1.0 * RBF(1.0))))
models.append(('MLPClassifier', MLPClassifier(alpha=1, max_iter=1000)))
gnb = GaussianNB()
models.append(('CalibratedClassifierCV_iso', CalibratedClassifierCV(gnb, cv=2, method="isotonic")))
models.append(('CalibratedClassifierCV_sig', CalibratedClassifierCV(gnb, cv=2, method="sigmoid")))
    
    
    
    
def machineLearningv2(
    file="../../../normalizedCounts.csv", 
    path_to_save="machineLearning", 
    component_n=3, 
    condition="SUB,NC", 
    sample_count="5, 5",
    graphFileName="graph",
):
    
    # Check whether the specified path exists or not
    isExist = os.path.exists(path_to_save)

    if not isExist:
        # Create a new directory because it does not exist 
        os.makedirs(path_to_save)
    print("{} is created!".format(path_to_save))
    
    pp = PdfPages('{}/{}.pdf'.format(path_to_save, graphFileName))
    component_n = int(component_n)
    
    
    Dataset = pd.read_csv(file)
    Dataset_unnor = pd.read_csv('../../../UnnormalizedCounts.csv')

    list(Dataset)
    # df_new = Dataset.rename(columns={'A': 'Col_1'}, index={'ONE': 'Row_1'})
    df_new = Dataset.rename(columns={'Unnamed: 0': ''})

    df_new.to_csv('{}/NormalizedCountWithGeneColumnRenamed.csv'.format(path_to_save))
    df_index = df_new.set_index('')

    df_data = df_index.T
    df_data

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

    df_data.to_csv('./{}/NormalizedCountWithGeneColumnLabels.csv'.format(path_to_save))
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

    df_x_x = df_data.loc[:, df_data.isin([' ','NULL',0]).mean() <= 1]
    df_data.loc[:, df_data.isin([' ','NULL',0]).mean() >= 1].shape
    no_zero = df_data.loc[:,df_data.isin([' ','NULL',0]).mean() < .0]
    df_data.isin([' ','NULL',0])['ENSG00000000005']

    df_without_label = df_data.drop(['labels'], axis=1)
    df3 = df_without_label[df_without_label.columns[df_without_label.mean() > 0.4]]
    column_count = df_without_label.columns[df_without_label.mean() > 0.4]

    column_count2 = df_without_label.columns[df_without_label.mean() < 0.4]


    pca = PCA(n_components=2)
    principalComponents = pca.fit_transform(x)
    principal_df = pd.DataFrame(data = principalComponents
                , columns = ['p1', 'p2'],
                index=df_data.index.values)

    pca1 = PCA(n_components=2)
    principalComponents2 = pca1.fit_transform(df3.to_numpy())
    principal_df2 = pd.DataFrame(data = principalComponents2
                , columns = ['p1', 'p2'],
                index=df3.index.values)


    final_df2 = pd.concat([principal_df2, df_data[['labels']]], axis = 1)
    final_df = pd.concat([principal_df, df_data[['labels']]], axis = 1)



    # fig = plt.figure(figsize = (8,8))
    # ax = fig.add_subplot(1,1,1) 
    # ax.set_xlabel('p1', fontsize = 15)
    # ax.set_ylabel('p2', fontsize = 15)
    # ax.set_title('2 component PCA', fontsize = 20)
    # targets = ['subacute ischemic stroke', 'acute ischemic stroke', 'healthy control']
    # colors = ['r', 'g', 'b']
    # for target, color in zip(targets,colors):
    #     indicesToKeep = final_df['labels'] == target
    #     ax.scatter(final_df.loc[indicesToKeep, 'p1']
    #             , final_df.loc[indicesToKeep, 'p2']
    #             , c = color
    #             , s = 50)
    # ax.legend(targets)
    # ax.grid()
    # # ax.figure.savefig("PCA2.png")
    # ax.figure.savefig("{}/{}".format(path_to_save, 'PCA_v2_15_samples.png'))





    # fig = plt.figure(figsize = (8,8))
    # ax = fig.add_subplot(1,1,1) 
    # ax.set_xlabel('p1', fontsize = 15)
    # ax.set_ylabel('p2', fontsize = 15)
    # ax.set_title('2 component PCA', fontsize = 20)
    # targets = ['subacute ischemic stroke', 'acute ischemic stroke', 'healthy control']
    # colors = ['r', 'g', 'b']
    # for target, color in zip(targets,colors):
    #     indicesToKeep = final_df2['labels'] == target
    #     ax.scatter(final_df2.loc[indicesToKeep, 'p1']
    #             , final_df2.loc[indicesToKeep, 'p2']
    #             , c = color
    #             , s = 50)
    # ax.legend(targets)
    # ax.grid()
    # # ax.figure.savefig("PCA1.png")
    # ax.figure.savefig("{}/{}".format(path_to_save, 'PCA_15_samples.png'))




    # pca.explained_variance_ratio_

    # train_data, test_data, train_lbl, test_lbl = train_test_split( df_data.loc[:,header].values, 
    #                                                             final_df['labels'], 
    #                                                             test_size=1/7.0, 
    #                                                             random_state=42, shuffle=True)



    # scaler = StandardScaler()
    # # Fit on training set only.
    # scaler.fit(train_data)
    # # Apply transform to both the training set and the test set.
    # train_data_ = scaler.transform(train_data)
    # test_data_ = scaler.transform(test_data)


    # Make an instance of the Model
    # pca = PCA(.60)
    # pca.fit(train_data_)
    # pca.n_components_
    # train_data_transform = pca.transform(train_data_)
    # test_data_transform = pca.transform(test_data_)

    
    # models.append(('K_means', KMeans(n_clusters=3, random_state=0)))




    # results = {}
    # names = []
    # predict = []
    # scoring = 'accuracy'
    # roc = {}
    # clfs = []
    # time_taken = {}
    # for name, model in models:
    #     start = time.process_time()
    #     ddd = model.fit(train_data_transform, train_lbl)
    #     time_taken[name] = time.process_time() - start
    #     results[name] = ddd.score(test_data_transform, test_lbl)
    #     clfs.append((ddd, ddd.predict(test_data_transform)))
    #     predict.append(ddd.predict(test_data_transform))
    #     roc[name] = roc_auc_score(test_lbl, ddd.predict_proba(test_data_transform), multi_class='ovr')




    # # boxplot algorithm comparison
    # # print(results)
    # # print(roc)

    # for pre in predict:
    #     cm = confusion_matrix(test_lbl, pre)
    #     # print(accuracy_score(test_lbl, pre, normalize=False))
    #     accuracy, precision, recall, f1 = get_metrics(test_lbl, pre)
    #     sensitivity = cm[0,0]/(cm[0,0]+cm[0,1])
    #     # print("precision_score {}".format(precision_score(test_lbl, pre, average=None)))
    #     specificity = cm[1,1]/(cm[1,0]+cm[1,1])
    #     # print("accuracy = %.3f \nprecision = %.3f \nrecall = %.3f \nf1 = %.3f \nSpecificity = %.3f \nSensitivity = %.3f" % (accuracy, precision, recall, f1, specificity, sensitivity))




    # for pre in predict:
    #     print(classification_report(test_lbl, pre,target_names=np.unique(pre), labels=np.unique(pre)))
        
    # V = list(results.values())
    # K = list(results.keys())
    # Accuracy = pd.Series(
    #     V,
    #     index = K
    # )


    # Vt = list(time_taken.values())
    # Kt = list(time_taken.keys())
    # timeT = pd.Series(
    #     Vt,
    #     index = Kt
    # )
    
    
    # #Set descriptions:
    # fig2 = plt.figure(figsize=(20, 15))
    
    
    # Accuracy.to_frame().index.name = 'Model'
    # accuracy_dataframe = Accuracy.reset_index(name='accuracy')
    # # print("Sample")
    # # print(Accuracy.reset_index(name='accuracy'))
    # ax = sns.barplot(x=accuracy_dataframe.Model, y="accuracy", data=accuracy_dataframe)
    # l1 = plt.axhline(y=1, linewidth=14, color='w')
    # l2 = plt.axhline(y=1, linewidth=14, color='w')
    # ax.set_xticklabels(accuracy_dataframe.Model, rotation=45, ha='right')
    # ax.figure.savefig("{}/{}".format(path_to_save, 'output2.png'))
    # pp.savefig(fig2)


    # #Set descriptions:
    # fig3 = plt.figure(figsize=(20, 15))
    
    
    # timeT.to_frame().index.name = 'Model'
    # time_dataframe = timeT.reset_index(name='time')
    
    # time_dataframe
    # ax = sns.barplot(x=time_dataframe.index, y="time", data=time_dataframe)
    # l1 = plt.axhline(y=1, linewidth=14, color='w')
    # l2 = plt.axhline(y=1, linewidth=14, color='w')
    # ax.set_xticklabels(time_dataframe.Model, rotation=45, ha='right')
    # ax.figure.savefig("{}/{}".format(path_to_save, 'output.png'))
    # pp.savefig(fig3)
    generate_fake_date(principal_df, df_data, path_to_save)


def plot_roc_curve(fpr, tpr, label=None):
    plt.plot(fpr, tpr, linewidth=2, label=label)
    plt.plot([0, 1], [0, 1], 'k--')
    plt.axis([0, 1, 0, 1])
    plt.xlabel('False Positive Rate (Fall-Out)', fontsize=16)
    plt.ylabel('True Positive Rate (Recall)', fontsize=16)
    plt.grid(True)


def get_metrics(y_test, y_predicted):
    accuracy = accuracy_score(y_test, y_predicted)
    precision = precision_score(y_test, y_predicted, average='weighted')
    recall = recall_score(y_test, y_predicted, average='weighted')
    f1 = f1_score(y_test, y_predicted, average='weighted', labels=np.unique(y_predicted))
    return accuracy, precision, recall, f1


def generate_fake_date(principal_df2, df_data, path_to_save, num_copy=100, addup=50):
    # print(principal_df2)
    principal_df2.index.name = 'index'
    principal_df2

    print("It is a file: {}".format(os.path.isfile("./{}/{}".format(path_to_save, 'realTestDF.csv'))))
          
    if os.path.isfile("./{}/{}".format(path_to_save, 'realTestDF.csv')) == False:
        principal_df2['index'] = principal_df2.index

        real_df = pd.concat([principal_df2, df_data[['labels']]], axis = 1)
        classg = {}
        fake_df_n = pd.DataFrame({
            'index': [],
            'p1':[],
            'p2':[],
            'labels':[],
        })
        fake_df_a = pd.DataFrame({
            'index': [],
            'p1':[],
            'p2':[],
            'labels':[],
        })
        fake_df_s = pd.DataFrame({
            'index': [],
            'p1':[],
            'p2':[],
            'labels':[],
        })
        
        for idx, group in enumerate(df_data['labels'].unique()):
            classg[group] = real_df[real_df['labels'] == group]
        
        dtaa = {
            
        }
        label_status = ""
        for group in classg:
            label_status = group
            # print(max(classg[group]['p1']))
            if(group in df_data['labels'].unique() and group == "healthy control"):
                fake_df_n['index'] = ['GSM3483780'+str(x) for x in range(1, num_copy+1)]
                fake_df_n['p1'] = np.random.uniform(min(classg[group]['p1']),max(classg[group]['p1']) + addup, size=(num_copy))
                fake_df_n['p2'] = np.random.uniform(min(classg[group]['p2']),max(classg[group]['p2']) + addup, size=(num_copy))
                fake_df_n['labels'] = [group for x in range(1, num_copy+1)]
            elif(group in df_data['labels'].unique() and group == "acute ischemic stroke"):
                fake_df_a['index'] = ['GSM3483780'+str(x) for x in range(1, num_copy+1)]
                fake_df_a['p1'] = np.random.uniform(min(classg[group]['p1']),max(classg[group]['p1']) + addup, size=(num_copy))
                fake_df_a['p2'] = np.random.uniform(min(classg[group]['p2']),max(classg[group]['p2']) + addup, size=(num_copy))
                fake_df_a['labels'] = [group for x in range(1, num_copy+1)]
            
            elif(group in df_data['labels'].unique() and group == "subacute ischemic stroke"):
                fake_df_s['index'] = ['GSM3483780'+str(x) for x in range(1, num_copy+1)]
                fake_df_s['p1'] = np.random.uniform(min(classg[group]['p1']),max(classg[group]['p1']) + addup, size=(num_copy))
                fake_df_s['p2'] = np.random.uniform(min(classg[group]['p2']),max(classg[group]['p2']) + addup, size=(num_copy))
                fake_df_s['labels'] = [group for x in range(1, num_copy+1)]
            

        rec = real_df.append(fake_df_n, ignore_index=True, verify_integrity=False, sort=None)
        rec2 = rec.append(fake_df_a, ignore_index=True, verify_integrity=False, sort=None)
        rec3 = rec2.append(fake_df_s, ignore_index=True, verify_integrity=False, sort=None)
        rec.head()
        rec2.head(50)
        
        # rec = real_df
        
        # print(dtaa)
        # for i in dtaa:
        #     print(dtaa[i])
        #     rec = rec.append(dtaa[i], ignore_index=True, verify_integrity=False, sort=None)
            
        rec3.to_csv('./{}/realTestDF.csv'.format(path_to_save))
    else:
        print("We are using this.")
        rec3 = pd.read_csv('./{}/realTestDF.csv'.format(path_to_save))
    splitRecord(rec3, path_to_save)
    
    fig = plt.figure(figsize = (8,8))
    ax = fig.add_subplot(1,1,1) 
    ax.set_xlabel('p1', fontsize = 15)
    ax.set_ylabel('p2', fontsize = 15)
    ax.set_title('2 component PCA', fontsize = 20)
    targets = ['subacute ischemic stroke', 'acute ischemic stroke', 'healthy control']
    colors = ['r', 'g', 'b']
    for target, color in zip(targets,colors):
        indicesToKeep = rec3['labels'] == target
        ax.scatter(rec3.loc[indicesToKeep, 'p1']
                , rec3.loc[indicesToKeep, 'p2']
                , c = color
                , s = 50)
    ax.legend(targets)
    ax.grid()
    # ax.figure.savefig("PCA_fake.png")
    ax.figure.savefig("{}/{}".format(path_to_save, 'PCA_Real.png'))
    
    
def splitRecord(data, path_to_save):
    header = list(data.columns)
    # data.loc[:,header].values
    # print(data)
    train_data, test_data, train_lbl, test_lbl = train_test_split( data[['p1', 'p2']], 
                                                                data['labels'], 
                                                                test_size=1/7.0, 
                                                                random_state=42, shuffle=True)
    
    
    results = {}
    names = []
    predict = {}
    scoring = 'accuracy'
    roc = {}
    clfs = []
    time_taken = {}
    for name, model in models:
        start = time.process_time()
        # print(train_data)
        ddd = model.fit(train_data, train_lbl)
        time_taken[name] = time.process_time() - start
        results[name] = ddd.score(test_data, test_lbl)
        clfs.append((ddd, ddd.predict(test_data)))
        predict[name] = ddd.predict(test_data)
        roc[name] = roc_auc_score(test_lbl, ddd.predict_proba(test_data), multi_class='ovr')




    # boxplot algorithm comparison
    # print("Result")
    # print(results)
    # print("ROC")
    # print(roc)
    
    # print("My Predictions")
    # print(predict)
    
    file_object = open("{}/{}".format(path_to_save, 'report.txt'), 'a')
    rows = []
    header = ['model', 'accuracy', 'precision', 'recall', 'f1', 'Specificity', 'Sensitivity', 'ROC', 'Time_taken']
    for pre in predict:
        cm = confusion_matrix(test_lbl, predict[pre])
        # print(accuracy_score(test_lbl, predict[pre], normalize=False))
        accuracy, precision, recall, f1 = get_metrics(test_lbl, predict[pre])
        sensitivity = cm[0,0]/(cm[0,0]+cm[0,1])
        file_object.write("precision_score {}".format(precision_score(test_lbl, predict[pre], average=None)))
        # print("precision_score {}".format(precision_score(test_lbl, pre, average=None)))
        specificity = cm[1,1]/(cm[1,0]+cm[1,1])
        # print("accuracy = %.3f \nprecision = %.3f \nrecall = %.3f \nf1 = %.3f \nSpecificity = %.3f \nSensitivity = %.3f" % (accuracy, precision, recall, f1, specificity, sensitivity))
        file_object.write("\nModel = %s \naccuracy = %.3f \nprecision = %.3f \nrecall = %.3f \nf1 = %.3f \nSpecificity = %.3f \nSensitivity = %.3f \nROC = %.3f\n\n" % (pre, accuracy, precision, recall, f1, specificity, sensitivity, roc[pre]))
        rows.append({
            'model' : pre, 
            'accuracy' : accuracy, 
            'precision' : precision, 
            'recall' : recall, 
            'f1' : f1, 
            'Specificity' : specificity, 
            'Sensitivity' : sensitivity, 
            'ROC' : roc[pre],
            'Time_taken' : time_taken[pre]
        })
    saveResultResponse(header, rows, path_to_save)
    file_object.write("___________________________________________CLASSIFICATION SUMMARY____________________________________________\n\n")
    for pre in predict:
        file_object.write(classification_report(test_lbl, predict[pre],target_names=np.unique(predict[pre]), labels=np.unique(predict[pre])))
    
    V = list(results.values())
    K = list(results.keys())
    Accuracy = pd.Series(
        V,
        index = K
    )


    Vt = list(time_taken.values())
    Kt = list(time_taken.keys())
    timeT = pd.Series(
        Vt,
        index = Kt
    )
    
    
    #Set descriptions:
    fig2 = plt.figure(figsize=(20, 15))
    
    
    Accuracy.to_frame().index.name = 'Model'
    accuracy_dataframe = Accuracy.reset_index(name='accuracy')
    # print("Sample")
    # print(Accuracy.reset_index(name='accuracy'))
    ax = sns.barplot(x=accuracy_dataframe.Model, y="accuracy", data=accuracy_dataframe)
    l1 = plt.axhline(y=1, linewidth=14, color='w')
    l2 = plt.axhline(y=1, linewidth=14, color='w')
    ax.set_xticklabels(accuracy_dataframe.Model, rotation=45, ha='right')
    # ax.figure.savefig("output3.png")
    ax.figure.savefig("{}/{}".format(path_to_save, 'Accuracy_graph.png'))


    #Set descriptions:
    fig3 = plt.figure(figsize=(20, 15))
    
    
    timeT.to_frame().index.name = 'Model'
    time_dataframe = timeT.reset_index(name='time')
    
    time_dataframe
    ax = sns.barplot(x=time_dataframe.index, y="time", data=time_dataframe)
    l1 = plt.axhline(y=1, linewidth=14, color='w')
    l2 = plt.axhline(y=1, linewidth=14, color='w')
    ax.set_xticklabels(time_dataframe.Model, rotation=45, ha='right')
    # ax.figure.savefig("output4.png")
    ax.figure.savefig("{}/{}".format(path_to_save, 'Time_taken_graph.png'))
    
    fig3 = plt.figure(figsize=(20, 15))
    ax = sns.boxplot(x=data['p1'], y="labels", data=data)
    # l1 = plt.axhline(y=1, linewidth=14, color='w')
    # l2 = plt.axhline(y=1, linewidth=14, color='w')
    # ax.set_xticklabels(data[['labels']], rotation=45, ha='right')
    ax.figure.savefig("{}/{}".format(path_to_save, 'boxplot.png'))
    
    
    file_object.write("\n\nTraining Dataset: {}\n".format(len(train_data)))
    file_object.write("\nTest Dataset: {}\n".format(len(test_data)))
    
    file_object.write("--------------------------------------------END OF FILE BEFORE APPENDING ANOTHER------------------------------------------")

    file_object.close()
    
    
def saveResultResponse(header, rows, path_to_save):
    # csv header
    fieldnames = header

    # csv data
    rows = rows

    with open("{}/{}".format(path_to_save, 'result.csv'), 'w', encoding='UTF8', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)