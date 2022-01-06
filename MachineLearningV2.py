import matplotlib.pyplot as plt
import pandas as pd
from sklearn.preprocessing import StandardScaler

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


Dataset = pd.read_csv('../../../normalizedCounts.csv')
Dataset_unnor = pd.read_csv('../../../UnnormalizedCounts.csv')

list(Dataset)
# df_new = Dataset.rename(columns={'A': 'Col_1'}, index={'ONE': 'Row_1'})
df_new = Dataset.rename(columns={'Unnamed: 0': ''})

df_new.to_csv('NormalizedCountWithGeneColumnRenamed.csv')
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

df_data.to_csv('NormalizedCountWithGeneColumnLabels.csv')
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

pca = PCA(n_components=2)
principalComponents2 = pca.fit_transform(df3.to_numpy())
principal_df2 = pd.DataFrame(data = principalComponents2
             , columns = ['p1', 'p2'],
            index=df3.index.values)


final_df2 = pd.concat([principal_df2, df_data[['labels']]], axis = 1)
final_df = pd.concat([principal_df, df_data[['labels']]], axis = 1)



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





fig = plt.figure(figsize = (8,8))
ax = fig.add_subplot(1,1,1) 
ax.set_xlabel('p1', fontsize = 15)
ax.set_ylabel('p2', fontsize = 15)
ax.set_title('2 component PCA', fontsize = 20)
targets = ['subacute ischemic stroke', 'acute ischemic stroke', 'healthy control']
colors = ['r', 'g', 'b']
for target, color in zip(targets,colors):
    indicesToKeep = final_df2['labels'] == target
    ax.scatter(final_df2.loc[indicesToKeep, 'p1']
               , final_df2.loc[indicesToKeep, 'p2']
               , c = color
               , s = 50)
ax.legend(targets)
ax.grid()




pca.explained_variance_ratio_


from sklearn.model_selection import train_test_split
train_data, test_data, train_lbl, test_lbl = train_test_split( df_data.loc[:,header].values, 
                                                              final_df['labels'], 
                                                              test_size=1/7.0, 
                                                              random_state=42, shuffle=True)



scaler = StandardScaler()
# Fit on training set only.
scaler.fit(train_data)
# Apply transform to both the training set and the test set.
train_data_ = scaler.transform(train_data)
test_data_ = scaler.transform(test_data)







from sklearn.decomposition import PCA
# Make an instance of the Model
pca = PCA(.60)
pca.fit(train_data_)
pca.n_components_
train_data_transform = pca.transform(train_data_)
test_data_transform = pca.transform(test_data_)







# logisticRegr = LogisticRegression(solver = 'lbfgs')
# pipe = make_pipeline(StandardScaler(), LogisticRegression())
# pipe_svm = make_pipeline(StandardScaler(), SVC())

# pipe_rdf = make_pipeline(StandardScaler(), RandomForestClassifier())
# pipe_lda = make_pipeline(StandardScaler(), LinearDiscriminantAnalysis())

# neigh = KNeighborsClassifier(n_neighbors=3)
# pipe_knn = make_pipeline(StandardScaler(), neigh)

# decision_tree = tree.DecisionTreeClassifier()
# pipe_dt = make_pipeline(StandardScaler(), decision_tree)

# kmeans = KMeans(n_clusters=2, random_state=0)
# pipe_km = make_pipeline(StandardScaler(), kmeans)


models = []
models.append(('LR', LogisticRegression(solver='lbfgs', max_iter=1000)))
models.append(('LDA', LinearDiscriminantAnalysis()))
models.append(('KNN', KNeighborsClassifier(n_neighbors=3)))
models.append(('CART', tree.DecisionTreeClassifier(max_depth=5)))
models.append(('NB', GaussianNB()))
models.append(('SVM', SVC(random_state=0, probability=True)))
models.append(('RandomForest', RandomForestClassifier(max_depth=5, n_estimators=10, max_features=1)))
models.append(('AdaBooster', AdaBoostClassifier(n_estimators=100, random_state=0)))
models.append(('BaggingClassifier', BaggingClassifier(base_estimator=SVC(), n_estimators=10, random_state=0)))
# models.append(('RidgeClassifierCV', RidgeClassifierCV(alphas=[1e-3, 1e-2, 1e-1, 1])))
models.append(('QuadraticDiscriminantAnalysis', QuadraticDiscriminantAnalysis()))
models.append(('SVC_linear', SVC(kernel="linear", C=0.025, probability=True)))
models.append(('SVC_gamma_2', SVC(gamma=2, C=1, probability=True)))
models.append(('GaussianProcessClassifier', GaussianProcessClassifier(1.0 * RBF(1.0))))
models.append(('MLPClassifier', MLPClassifier(alpha=1, max_iter=1000)))
gnb = GaussianNB()
models.append(('CalibratedClassifierCV_iso', CalibratedClassifierCV(gnb, cv=2, method="isotonic")))
models.append(('CalibratedClassifierCV_sig', CalibratedClassifierCV(gnb, cv=2, method="sigmoid")))
# models.append(('K_means', KMeans(n_clusters=3, random_state=0)))




results = {}
names = []
predict = []
scoring = 'accuracy'
roc = {}
clfs = []
for name, model in models:
    ddd = model.fit(train_data_transform, train_lbl)
    results[name] = ddd.score(test_data_transform, test_lbl)
    clfs.append((ddd, ddd.predict(test_data_transform)))
    predict.append(ddd.predict(test_data_transform))
    roc[name] = roc_auc_score(test_lbl, ddd.predict_proba(test_data_transform), multi_class='ovr')
#     kfold = model_selection.KFold(n_splits=1, random_state=seed)
#     cv_results = model_selection.cross_val_score(model, test_data_transform, test_lbl, cv=kfold, scoring=scoring)
#     cv_results = model.fit(test_data_transform, test_lbl)
#     results.append(cv_results)
#     names.append(name)
#     print(results)
#     msg = "%s: %f (%f)" % (name, cv_results.mean(), cv_results.std())
#     print(msg)




# boxplot algorithm comparison
print(results)
print(roc)

# fig = plt.figure()
# fig.suptitle('Algorithm Comparison')
# ax = fig.add_subplot(111)
# plt.boxplot(results)
# ax.set_xticklabels(names)
# plt.show()
from sklearn.metrics import accuracy_score, confusion_matrix, f1_score, precision_score, recall_score
from sklearn import metrics
import numpy as np

def get_metrics(y_test, y_predicted):
    accuracy = accuracy_score(y_test, y_predicted)
    precision = precision_score(y_test, y_predicted, average='weighted')
    recall = recall_score(y_test, y_predicted, average='weighted')
    f1 = f1_score(y_test, y_predicted, average='weighted', labels=np.unique(y_predicted))
    return accuracy, precision, recall, f1
for pre in predict:
    cm = confusion_matrix(test_lbl, pre)
    print(accuracy_score(test_lbl, pre, normalize=False))
    accuracy, precision, recall, f1 = get_metrics(test_lbl, pre)
    sensitivity = cm[0,0]/(cm[0,0]+cm[0,1])
    print("precision_score {}".format(precision_score(test_lbl, pre, average=None)))
    specificity = cm[1,1]/(cm[1,0]+cm[1,1])
    print("accuracy = %.3f \nprecision = %.3f \nrecall = %.3f \nf1 = %.3f \nSpecificity = %.3f \nSensitivity = %.3f" % (accuracy, precision, recall, f1, specificity, sensitivity))




for pre in predict:
    print(classification_report(test_lbl, pre,target_names=np.unique(pre), labels=np.unique(pre)))
