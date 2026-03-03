#!/usr/bin/env python3

"""
Fitting and testing various models for binary classification

INPUTS
------

data
    Tab-separated table of shape N * (M + 1), where N is the number of samples and M are the
    variables to be used for prediction.
    The first column should be unique sample identifiers.
    The first row should be a header.

labels
    Tab-delimited file with two columns: samples IDs (matching the data table) and labels for
    predictions.
    Expected to NOT have a header row.
"""

### ---------------------------------------- ###

def parse_args():

    ### Data
    data_file = argv[argv.index('--data') + 1]
    data = pd.read_csv(data_file, sep='\t', index_col=0, header=0)

    ### Labels
    labels_file = argv[argv.index('--labels') + 1]
    labels = pd.read_csv(labels_file, sep='\t', index_col=0, header=None)
    
    ### Sanity check
    good_samples = data.index.values[data.index.isin(labels.index.values)]
    data = data.loc[good_samples,]
    labels = labels.loc[good_samples,]

    ### Rename labels
    unique_labels = {l : i for i,l in enumerate(np.unique(labels.values))}
    labels = labels.replace(unique_labels)

    return data, labels, unique_labels

### ---------------------------------------- ###

def test_model(X, y, regr, name, split_call):
    
    print('#' * 40)
    
    print(f'### Evaluating {name}')
    
    # Init results collection
    model_results = {'model' : name,
                     'train_mse_mean' : [],
                     'train_mse_std' : 0,
                     'train_r2_mean' : [],
                     'train_r2_std' : 0,
                     'train_overall_accuracy_mean' : [],
                     'train_overall_accuracy_std' : 0,
                     'train_class-0_accuracy_mean' : [],
                     'train_class-0_accuracy_std' : 0,
                     'train_class-1_accuracy_mean' : [],
                     'train_class-1_accuracy_std' : 0,
                     'test_mse_mean' : [],
                     'test_mse_std' : 0,
                     'test_r2_mean' : [],
                     'test_r2_std' : 0,
                     'test_overall_accuracy_mean' : [],
                     'test_overall_accuracy_std' : 0,
                     'test_class-0_accuracy_mean' : [],
                     'test_class-0_accuracy_std' : 0,
                     'test_class-1_accuracy_mean' : [],
                     'test_class-1_accuracy_std' : 0}
    
    # Cross validation
    for n, (train_indexes, test_indexes) in enumerate(split_call.split(X, y)):
        
        print(f'## Training {n + 1}:')
        
        # Split data
        train_x, train_y = X.iloc[train_indexes].copy(), y.values.flatten()[train_indexes].copy()
        test_x, test_y = X.iloc[test_indexes].copy(), y.values.flatten()[test_indexes].copy()
        
        # Copy model without fit
        regr_test = clone(regr)

        # Fit data
        regr_test.fit(train_x, train_y)
        
        # Test model
        for dataset_name,x_dat,y_dat in zip(['train', 'test'], [train_x, test_x], [train_y, test_y]):
            
            # Assess model
            prediction = regr_test.predict(x_dat)
            mse = mean_squared_error(y_dat, prediction)
            r2 = r2_score(y_dat, prediction)
            overall_correct = (y_dat == prediction).sum() / len(y_dat)
            category_0_tot = (y_dat == 0).sum()
            category_0_correct = ((y_dat == prediction) & (y_dat == 0)).sum() / category_0_tot
            category_1_tot = (y_dat == 1).sum()
            category_1_correct = ((y_dat == prediction) & (y_dat == 1)).sum() / category_1_tot
            
            # Store data
            model_results[f'{dataset_name}_mse_mean'].append(mse)
            model_results[f'{dataset_name}_r2_mean'].append(r2)
            model_results[f'{dataset_name}_overall_accuracy_mean'].append(overall_correct)
            model_results[f'{dataset_name}_class-0_accuracy_mean'].append(category_0_correct)
            model_results[f'{dataset_name}_class-1_accuracy_mean'].append(category_1_correct)
            
            print(f'# Accuracy {dataset_name} set:')
            print(f'Overall\t{(overall_correct * 100):.2F}%')
            print(f'Class 0\t{(category_0_correct * 100):.2F}%')
            print(f'Class 1\t{(category_1_correct * 100):.2F}%')
        
        print()
            
    # Average results and calculate standard deviations
    for dataset_name in ['train', 'test']:
        
        for param in ['mse', 'r2', 'overall_accuracy', 'class-0_accuracy', 'class-1_accuracy']:
            
            mean, std = np.mean(model_results[f'{dataset_name}_{param}_mean']), np.std(model_results[f'{dataset_name}_{param}_mean'])
            
            model_results[f'{dataset_name}_{param}_mean'] = mean
            model_results[f'{dataset_name}_{param}_std'] = std
    
    model_results = pd.DataFrame(model_results, index=[0])
    
    print('#' * 40)
    
    return model_results
    
### ------------------MAIN------------------ ###

import numpy as np
import pandas as pd
import seaborn as sns

from matplotlib import pyplot as plt
from sklearn.base import clone
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis, QuadraticDiscriminantAnalysis
from sklearn.ensemble import AdaBoostClassifier, RandomForestClassifier
from sklearn.linear_model import Lasso, LogisticRegression, Ridge
from sklearn.metrics import mean_squared_error, r2_score
from sklearn.model_selection import StratifiedShuffleSplit
from sklearn.naive_bayes import GaussianNB
from sklearn.neighbors import KNeighborsClassifier
from sklearn.neural_network import MLPClassifier
from sklearn.svm import SVC
from sys import argv

### Parse args

data, labels, unique_labels = parse_args()

### Init split between train and testing

# Using StratifiedShuffleSplit instead of RepeatedStratifiedKFold to account for imbalanced labels
sss = StratifiedShuffleSplit(n_splits=10, test_size=0.2, random_state=42)

### Model

models = {'AdaBoost' : AdaBoostClassifier(algorithm="SAMME", random_state=42),
          'GaussianNB' : GaussianNB(),
          'LinearDiscriminant' : LinearDiscriminantAnalysis(store_covariance=True),
          'QuadraticDiscriminant' : QuadraticDiscriminantAnalysis(store_covariance=True),
          'LogisticRegression' : LogisticRegression(fit_intercept=False, max_iter=10000),
          'Lasso' : Lasso(alpha=1, fit_intercept=False),
          'Ridge' : Ridge(alpha=1, fit_intercept=False),
          'NeuralNet' : MLPClassifier(alpha=1, max_iter=1000, random_state=42),
          'KNeighbors' : KNeighborsClassifier(3),
          'SVM' : SVC(kernel="linear", C=0.025, random_state=42),
          'RandomForest' : RandomForestClassifier(max_depth=5, n_estimators=10, max_features=1, random_state=42)}

models_evaluation = []

for model_name,model_call in models.items():
    
    res = test_model(data.copy(), labels.copy(), model_call, model_name, sss)
    
    models_evaluation.append(res)

# Concat results and save to file
models_evaluation = pd.concat(models_evaluation, axis=0, ignore_index=True)
models_evaluation.to_csv('ml_test.tsv', sep='\t', index=False)

### Plot data

# Prepare data for plotting
#models_evaluation = pd.read_csv('ml_test.tsv', sep='\t')
plot_data = models_evaluation.loc[:, ['model',
                                      'test_overall_accuracy_mean', 'test_overall_accuracy_std',
                                      'test_class-0_accuracy_mean', 'test_class-0_accuracy_std',
                                      'test_class-1_accuracy_mean', 'test_class-1_accuracy_std']].copy()
plot_data.columns = [c.replace('test_', '').replace('_', '\n') for c in plot_data.columns]

plot_data = pd.DataFrame({'x' : np.concatenate([plot_data.model for _ in range(3)]),
                          'y' : np.concatenate([plot_data[f'{acc_type}\naccuracy\nmean'].values for acc_type in ['overall', 'class-0', 'class-1']]),
                          'y_std' : np.concatenate([plot_data[f'{acc_type}\naccuracy\nstd'].values for acc_type in ['overall', 'class-0', 'class-1']]),
                          'Accuracy' : [acc_type for acc_type in ['Overall', 'class-0', 'class-1'] for _ in range(plot_data.shape[0])]})

# Plot
plt.figure(figsize=(10, 6))
sns.barplot(data=plot_data,
            x='x',
            y='y',
            hue='Accuracy',
            errorbar=None)

for i in range(int(plot_data.shape[0] / 3)): # Draw barplot
    
    means = models_evaluation.iloc[i,][[f'test_{acc_type}_accuracy_mean' for acc_type in ['overall', 'class-0', 'class-1']]].values
    stds = models_evaluation.iloc[i,][[f'test_{acc_type}_accuracy_std' for acc_type in ['overall', 'class-0', 'class-1']]].values
    mins, maxs = means - stds, means + stds
    xs = [i - 0.25, i, i + 0.25]
    plt.vlines(xs, mins, maxs, colors='black')

plt.xticks(rotation=90)
legend = plt.legend(bbox_to_anchor=(1, 1), loc='best')
plt.tight_layout()
plt.savefig('ml_accuracy_test.png', dpi=300)
plt.close()
