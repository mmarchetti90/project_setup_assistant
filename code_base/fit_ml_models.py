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

N.B. Use the "--binary" toggle for binary classification problems
"""

### IMPORTS -------------------------------- ###

import numpy as np
import pandas as pd
import seaborn as sns

from collections.abc import Callable
from datetime import datetime
from matplotlib import pyplot as plt
from sklearn.base import clone
from sklearn.discriminant_analysis import (
    LinearDiscriminantAnalysis,
    QuadraticDiscriminantAnalysis
)
from sklearn.ensemble import (
    AdaBoostClassifier,
    AdaBoostRegressor,
    BaggingClassifier,
    BaggingRegressor,
    ExtraTreesClassifier,
    ExtraTreesRegressor,
    GradientBoostingClassifier,
    GradientBoostingRegressor,
    HistGradientBoostingClassifier,
    HistGradientBoostingRegressor,
    RandomForestClassifier,
    RandomForestRegressor
)
from sklearn.gaussian_process import (
    GaussianProcessClassifier,
    GaussianProcessRegressor
)
from sklearn.linear_model import (
    ARDRegression,
    BayesianRidge,
    ElasticNet,
    GammaRegressor,
    HuberRegressor,
    Lars,
    Lasso,
    LassoLars,
    LinearRegression,
    LogisticRegression,
    MultiTaskElasticNet,
    MultiTaskLasso,
    OrthogonalMatchingPursuit,
    PoissonRegressor,
    QuantileRegressor,
    RANSACRegressor,
    Ridge,
    SGDClassifier,
    SGDRegressor,
    TheilSenRegressor,
    TweedieRegressor
)
from sklearn.metrics import (
    mean_squared_error,
    r2_score
)
from sklearn.model_selection import (
    KFold,
    StratifiedShuffleSplit
)
from sklearn.naive_bayes import (
    BernoulliNB,
    GaussianNB,
)
from sklearn.neighbors import KNeighborsClassifier
from sklearn.neural_network import (
    MLPClassifier,
    MLPRegressor
)
from sklearn.svm import (
    SVC,
    SVR
)
from sys import argv

### CLASSES -------------------------------- ###

class model_test:
    
    """
    Class for evaluating model
    
    model_name : str
        name of the model
    model : Callable
        Model to be evaluated
    
    Methods
    -------
    forward(self, data: pd.DataFrame, split: Callable, binary: bool=False) -> dict[str, str, pd.DataFrame]
        Fits and evaluates the model
    log(log_message: str, timestamp: bool=True) -> None
        Prints a log message to stdout
    """
    
    def __init__(self, model_name: str, model: Callable) -> None:
        
        """
        Class init
        
        Parameters
        ----------
        model_name : str
            name of the model
        model : Callable
            Model to be evaluated
        """
        
        self.model_name = model_name
        
        self.model = model

    ### ------------------------------------ ###
    
    @staticmethod
    def log(log_message: str, timestamp: bool=True) -> None:
        
        """
        Prints a log message to stdout
        
        Parameters
        ----------
        log_message : str
            Message to print
        timestamp : bool
            Set to True to prepend a timestamp
        """
        
        if timestamp:
    
            log_timestamp = datetime.now().strftime("%H:%M:%S")
    
            timestamped_log_message = f"[{log_timestamp}] {log_message}"
    
            print(timestamped_log_message)
        
        else:
            
            print(log_message)
    
    ### ------------------------------------ ###
    
    def forward(self, X: pd.DataFrame, y: pd.DataFrame, split_call: Callable, binary: bool=False) -> pd.DataFrame:
        
        """
        Fits and evaluates the model
        
        Parameters
        ----------
        X : pd.DataFrame
            Variables to be used for prediction
        y : pd.DataFrame
            Variable to predict
        split_call : Callable
            K-fold splitting call
        binary : bool=False
            Set to True to evaluate binary classification
            Set to False for continuous variables
    
        Returns
        -------
        pd.DataFrame
            DataFrame model evaluation report
        """
        
        self.log(f'Evaluating {self.model_name}')
        
        # Init results collection
        
        if not binary:
            
            model_results = {
                'model' : [],
                'test_n' : [],
                'train_mse' : [],
                'train_r2' : [],
                'test_mse' : [],
                'test_r2' : []
            }
        
        else:
        
            model_results = {
                'model' : [],
                'test_n' : [],
                'train_mse' : [],
                'train_r2' : [],
                'train_overall_accuracy' : [],
                'train_class-0_accuracy' : [],
                'train_class-1_accuracy' : [],
                'test_mse' : [],
                'test_r2' : [],
                'test_overall_accuracy' : [],
                'test_class-0_accuracy' : [],
                'test_class-1_accuracy' : [],
            }
        
        # Cross validation
        
        for n, (train_indexes, test_indexes) in enumerate(split_call.split(X, y)):
            
            self.log(f'  * Training test {n + 1}')

            successes = 0
            
            try:
            
                # Split data into training and validation sets
                
                train_x, train_y = X.iloc[train_indexes].copy(), y.values.flatten()[train_indexes].copy()
                test_x, test_y = X.iloc[test_indexes].copy(), y.values.flatten()[test_indexes].copy()
                
                # Copy model without fit
    
                model_test = clone(self.model)
    
                # Fit data
                
                model_test.fit(train_x, train_y)

                model_results['model'].append(self.model_name)
                
                model_results['test_n'].append(n)
                
                successes += 1
                
                self.log('  \u2714 Training complete')
                
                # Test model
                
                for dataset_name,x_dat,y_dat in zip(['train', 'test'], [train_x, test_x], [train_y, test_y]):
                    
                    # Assess model
                    
                    prediction = model_test.predict(x_dat)
                    mse = mean_squared_error(y_dat, prediction)
                    r2 = r2_score(y_dat, prediction)
                    
                    # Store data
                    
                    model_results[f'{dataset_name}_mse'].append(mse)
                    model_results[f'{dataset_name}_r2'].append(r2)
                    
                    # Additional parameters for binary classification
                    
                    if binary:
                        
                        overall_correct = (y_dat == prediction).sum() / len(y_dat)
                        category_0_tot = (y_dat == 0).sum()
                        category_0_correct = ((y_dat == prediction) & (y_dat == 0)).sum() / category_0_tot
                        category_1_tot = (y_dat == 1).sum()
                        category_1_correct = ((y_dat == prediction) & (y_dat == 1)).sum() / category_1_tot
                        
                        # Store data
                        
                        model_results[f'{dataset_name}_overall_accuracy'].append(overall_correct)
                        model_results[f'{dataset_name}_class-0_accuracy'].append(category_0_correct)
                        model_results[f'{dataset_name}_class-1_accuracy'].append(category_1_correct)
                        
                        self.log(f'    * Accuracy on {dataset_name} set:')
                        self.log(f'      Overall\t{(overall_correct * 100):.2F}%')
                        self.log(f'      Class 0\t{(category_0_correct * 100):.2F}%')
                        self.log(f'      Class 1\t{(category_1_correct * 100):.2F}%')
                
            except:
                    
                self.log('  \u2717 Training failed')
        
        model_results = pd.DataFrame(model_results)
        
        # Complete message
        
        self.log('  \u2714 Evaluation complete')
        
        return model_results, successes > 0

### FUNCTIONS ------------------------------ ###

def parse_args():
    
    """
    Parsing CLI args
    """

    # Data
    
    data_file = argv[argv.index('--data') + 1]
    data = pd.read_csv(data_file, sep='\t', index_col=0, header=0)

    # Labels
    
    labels_file = argv[argv.index('--labels') + 1]
    labels = pd.read_csv(labels_file, sep='\t', index_col=0, header=None)
    
    # Mode
    
    binary_mode = ('--binary' in argv)
    
    # Sanity check
    
    good_samples = data.index.values[data.index.isin(labels.index.values)]
    data = data.loc[good_samples,]
    labels = labels.loc[good_samples,]

    # Rename labels
    
    unique_labels = {l : i for i,l in enumerate(np.unique(labels.values))}
    labels = labels.replace(unique_labels)

    return data, labels, unique_labels, binary_mode

### ---------------------------------------- ###

def init_models(binary: bool=False) -> tuple[Callable]:
    
    """
    Init a set of models to be tested
    
    Parameters
    ----------
    binary:
        Set to True to evaluate binary classification
        Set to False for continuous variables
    
    Returns
    -------
    models : dict[Callable]
        Tuple of model_test classes
    """
    
    if binary:
        
        models = (
            model_test('LogisticRegression', LogisticRegression(fit_intercept=False, max_iter=10000)),
            model_test('SGD', SGDClassifier()),
            model_test('NeuralNet', MLPClassifier(alpha=1, max_iter=1000, random_state=42)),
            model_test('LinearDiscriminant', LinearDiscriminantAnalysis(store_covariance=True)),
            model_test('QuadraticDiscriminant', QuadraticDiscriminantAnalysis(store_covariance=True)),
            model_test('AdaBoost', AdaBoostClassifier(algorithm="SAMME", random_state=42)),
            model_test('Bagging', BaggingClassifier(random_state=42)),
            model_test('ExtraTrees', ExtraTreesClassifier(random_state=42)),
            model_test('GradientBoosting', GradientBoostingClassifier(random_state=42)),
            model_test('HistGradientBoosting', HistGradientBoostingClassifier(random_state=42)),
            model_test('RandomForest', RandomForestClassifier(max_depth=5, n_estimators=10, max_features=1, random_state=42)),
            model_test('GaussianProcess', GaussianProcessClassifier(random_state=42)),
            model_test('BernoulliNB', BernoulliNB()),
            model_test('GaussianNB', GaussianNB()),
            model_test('KNeighbors', KNeighborsClassifier(3)),
            model_test('SVM', SVC(kernel="linear", C=0.025, random_state=42))
        )
    
    else:
    
        models = (
            model_test('ARD', ARDRegression()),
            model_test('BayesianRidge', BayesianRidge()),
            model_test('ElasticNet', ElasticNet(alpha=1, fit_intercept=False)),
            model_test('GammaRegressor', GammaRegressor(alpha=1, fit_intercept=False)),
            model_test('HuberRegressor', HuberRegressor(alpha=1, fit_intercept=False)),
            model_test('Lars', Lars(fit_intercept=False)),
            model_test('Lasso', Lasso(alpha=1, fit_intercept=False)),
            model_test('LassoLars', LassoLars(alpha=1, fit_intercept=False)),
            model_test('LinearRegression', LinearRegression(fit_intercept=False)),
            model_test('MultiTaskElasticNet', MultiTaskElasticNet(alpha=1, fit_intercept=False)),
            model_test('MultiTaskLasso', MultiTaskLasso(alpha=1, fit_intercept=False)),
            model_test('OMP', OrthogonalMatchingPursuit(fit_intercept=False)),
            model_test('Poisson', PoissonRegressor(alpha=1, fit_intercept=False)),
            model_test('Quantile', QuantileRegressor(alpha=1, fit_intercept=False)),
            model_test('RANSAC', RANSACRegressor()),
            model_test('Ridge', Ridge(alpha=1, fit_intercept=False)),
            model_test('SGD', SGDRegressor(fit_intercept=False)),
            model_test('TheilSen', TheilSenRegressor(fit_intercept=False)),
            model_test('Tweedie', TweedieRegressor(alpha=1, fit_intercept=False)),
            model_test('AdaBoost', AdaBoostRegressor(random_state=42)),
            model_test('Bagging', BaggingRegressor(random_state=42)),
            model_test('ExtraTrees', ExtraTreesRegressor(random_state=42)),
            model_test('GradientBoosting', GradientBoostingRegressor(random_state=42)),
            model_test('HistGradientBoosting', HistGradientBoostingRegressor(random_state=42)),
            model_test('GaussianProcess', GaussianProcessRegressor(random_state=42)),
            model_test('RandomForest', RandomForestRegressor(max_depth=5, n_estimators=10, max_features=1, random_state=42)),
            model_test('NeuralNet', MLPRegressor(alpha=1, max_iter=1000, random_state=42)),
            model_test('SVM', SVR(kernel="linear", C=0.025))
        )
    
    return models

### ---------------------------------------- ###

def evaluate_all_models(X: pd.DataFrame, y: pd.DataFrame, models: tuple[Callable], binary: bool=False):
    
    """
    Wrapper for sequentially testing each model and compiling a report
    
    Parameters
    ----------
    X : pd.DataFrame
        Variables to be used for prediction
    y : pd.DataFrame
        Variable to predict
    models : set[Callable]
        Tuple of model_test classes
    binary:
        Set to True to evaluate binary classification
        Set to False for continuous variables
    """
    
    # Init split between train and testing

    if binary:

        # Using StratifiedShuffleSplit instead of RepeatedStratifiedKFold to account for imbalanced labels
        
        split = StratifiedShuffleSplit(n_splits=10, test_size=0.2, random_state=42)

    else:
        
        split = KFold(n_splits=10, shuffle=True, random_state=42)
        
    # Evaluate models
    
    models_evaluation = []

    for model_class in models:
        
        res, ok = model_class.forward(X.copy(), y.copy(), split, binary)
        
        if ok:
        
            models_evaluation.append(res)

    # Concat results and save to file
    
    models_evaluation = pd.concat(models_evaluation, axis=0, ignore_index=True)
    
    return models_evaluation

### ---------------------------------------- ###

def plot_continuous_evaluation(models_evaluation: pd.DataFrame, output_name: str='ml_evaluation_test.png'):
    
    """
    Plots MSE and R2 for continuous predictions
    
    Parameters
    ----------
    models_evaluation : pd.DataFrame
        Evaluation results
    output_name : str
        Name of plot file
    """
    
    # Prepare data for plotting
    
    plot_data = models_evaluation.loc[:, ['model',
                                          'test_mse',
                                          'test_r2']].copy()
    
    plot_data.columns = [c.replace('test_', '').replace('_', '\n') for c in plot_data.columns]
    
    # Plot
    
    fig, axes = plt.subplots(nrows=2, ncols=1, figsize=(np.unique(plot_data['model']).shape[0], 12), sharex=True)
    
    for ax,param in zip(axes, ['mse', 'r2']):
        
        # Set y limits
        
        min_y, max_y = plot_data[param].quantile(0.05), plot_data[param].quantile(0.95)
    
        # Boxplot
        
        sns.boxplot(
            data=plot_data,
            x='model',
            y=param,
            hue='model',
            ax=ax
        )
        
        # Stripplot
        
        sns.stripplot(
            data=plot_data,
            x='model',
            y=param,
            color='black',
            ax=ax
        )
        
        ax.set_ylabel(param, fontweight='bold')
        
        ax.set_ylim(min_y, max_y)
    
    plt.xlabel(None)
    plt.xticks(rotation=45, fontweight='bold', ha='right')
    plt.tight_layout()
    plt.savefig(output_name, dpi=300)
    plt.close()

### ---------------------------------------- ###

def plot_binary_evaluation(models_evaluation: pd.DataFrame, output_name: str='ml_accuracy_test.png'):
    
    """
    Plots accuracy for binary classification
    
    Parameters
    ----------
    models_evaluation : pd.DataFrame
        Evaluation results
    output_name : str
        Name of plot file
    """
    
    # Prepare data for plotting
    
    plot_data = models_evaluation.loc[:, ['model',
                                          'test_overall_accuracy',
                                          'test_class-0_accuracy',
                                          'test_class-1_accuracy']].copy()
    
    plot_data.columns = [c.replace('test_', '').replace('_', '\n') for c in plot_data.columns]
    
    # Plot
    
    fig, axes = plt.subplots(nrows=3, ncols=1, figsize=(np.unique(plot_data['model']).shape[0], 12), sharex=True)
    
    for ax,param in zip(axes, ['overall\naccuracy', 'class-0\naccuracy', 'class-1\naccuracy']):
    
        sns.boxplot(
            data=plot_data,
            x='model',
            y=param,
            hue='model',
            ax=ax
        )
        
        sns.stripplot(
            data=plot_data,
            x='model',
            y=param,
            color='black',
            ax=ax
        )
        
        ax.set_ylabel(param, fontweight='bold')
    
    plt.xlabel(None)
    plt.xticks(rotation=45, fontweight='bold', ha='right')
    plt.tight_layout()
    plt.savefig(output_name, dpi=300)
    plt.close()

### MAIN ----------------------------------- ###

if __name__ == '__main__':

    # Parse args
    
    data, labels, unique_labels, binary_mode = parse_args()
    
    output_suffix = 'continuous' if not binary_mode else 'binary'
    
    # Init models set
    
    models = init_models(binary_mode)
    
    # Evaluate models
    
    models_evaluation = evaluate_all_models(data, labels, models, binary_mode)
    
    models_evaluation.to_csv(f'ml_test_{output_suffix}.tsv', sep='\t', index=False)
    
    # Plot data
    
    if not binary_mode:
        
        plot_continuous_evaluation(models_evaluation, 'ml_continuous_test')
    
    else:
    
        plot_binary_evaluation(models_evaluation, 'ml_binary_test.png')
