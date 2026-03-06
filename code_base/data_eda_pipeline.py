#!/usr/bin/env python3

"""
Exploratory Data Analysis pipeline
"""

### IMPORTS -------------------------------- ###

import numpy as np
import pandas as pd

from collections.abc import Callable
from datetime import datetime
from scipy.spatial import distance
from scipy.stats import (
    boxcox,
    chi2,
    chisquare,
    mannwhitneyu,
    normaltest
)
from warnings import catch_warnings, simplefilter

### CLASSES AND FUNCTIONS ------------------ ###

def help_fun():
    
    print("""
    Exploratory Data Analysis pipeline
    
    USAGE
    -----
        python data_eda_pipeline.py [args]
    
    ARGS
    ----
        --data_path
            Path to data file
            Required
        --sheet
            Sheet of xlsx file to be used
            Omit if file is text
        --onehot
            Include to one-hot encode categorical variables as a first step
        --remove_outliers
            Include to remove outliers
        --apply_transforms
            Include to apply transformations for improving normality
""")

### ---------------------------------------- ###

def parse_args():
    
    """
    Parses CLI arguments
    
    Returns
    -------
    data_path : str
        Path to data file
    sheet : str
        Sheet of xlsx file to be used
    onehot : bool
        If True, one-hot encodes categorical variables as a first step
    remove_outliers : bool
        If True, removes outlier samples
    apply_transforms : bool
        If True, applies transformations for improving normality
    """
    
    # Path to data (tsv, csv, txt, or xlsx table)
    
    data_path = argv[argv.index('--data_path') + 1] if '--data_path' in argv else ''
    
    # Sheet for xlsx file
    
    sheet = argv[argv.index('--sheet') + 1] if '--sheet' in argv else '0'
    
    # One-hot encoding
    
    one_hot = ('--onehot' in argv)
    
    # Remove outlier samples
    
    remove_outliers = ('--remove_outliers' in argv)
    
    # Apply transformations
    
    apply_transform = ('--apply_transforms' in argv)
    
    return data_path, sheet, one_hot, remove_outliers, apply_transform

### ---------------------------------------- ###

class pipeline:

    """
    Main class for analysis pipeline
    
    Parameters
    ----------
    pipeline_name : str
        Name used by the pipeline
    pipeline_description : str
        Concise description of the pipeline
    analysis_fun : Callable
        Function to apply to the dataset
    
    Methods
    -------
    forward(data: pd.DataFrame)
        Applies the analysis_fun to data
    """

    def __init__(self, pipeline_name: str, pipeline_description: str, analysis_fun: Callable) -> None:
        
        """
        Class init
        
        Parameters
        ----------
        pipeline_name : str
            Name used by the pipeline
        pipeline_description : str
            Concise description of the pipeline
        analysis_fun : Callable
            Function to apply to the dataset
        """

        self.pipeline_name = pipeline_name

        self.pipeline_description = pipeline_description

        self.analysis_fun = analysis_fun

    ### ------------------------------------ ###

    def forward(self, data: pd.DataFrame) -> tuple[pd.DataFrame, list[str]]:
        
        """
        Applies the analysis_fun to data
        
        Parameters
        ----------
        data : pd.DataFrame
            Dataframe to process
        
        Returns
        -------
        transformed_data : pd.DataFrame
            Data transformed after applying analysis_fun
        analysis_log : list[str]
            List of log messages
        """

        transformed_data, analysis_log = self.analysis_fun(data)

        return transformed_data, analysis_log

### ---------------------------------------- ###

class pipelines_manager:
    
    """
    Class for handling pipelines execution
    N.B. Pipelines are executed sequentially
    
    Parameters
    ----------
    pipeline_name : str
        Name used by the pipeline
    pipeline_description : str
        Concise description of the pipeline
    analysis_fun : Callable
        Function to apply to the dataset
    
    Methods
    -------
    log(log_message: str, timestamp: bool=True)
        Prints a log message to stdout
    add_pipeline(new_pipeline: Callable)
        Adds a pipeline class object
    explore(data: pd.DataFrame)
        Runs pipeline sequentially on the data
    """

    def __init__(self) -> None:
        
        """
        Class init
        
        Parameters
        ----------
        pipeline_name : str
            Name used by the pipeline
        pipeline_description : str
            Concise description of the pipeline
        analysis_fun : Callable
            Function to apply to the dataset
        """

        self.pipelines = []

        self.log(log_message='EDA pipelines manager initialized', timestamp=True)

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

    def add_pipeline(self, new_pipeline: Callable) -> None:
        
        """
        Adds a pipeline class object
        
        Parameters
        ----------
        new_pipeline : Callable
            New pipeline class object to add
        """

        self.pipelines.append(new_pipeline)

        self.log(log_message=f'Added new pipeline: {new_pipeline.pipeline_name}', timestamp=True)

    ### ------------------------------------ ###

    def explore(self, data: pd.DataFrame) -> dict[pd.DataFrame]:
        
        """
        Runs pipeline sequentially on the data
        
        Parameters
        ----------
        data : pd.DataFrame
            Dataframe to process
        
        Returns
        -------
        transformed_data : pd.DataFrame
            Data transformed after applying analysis_fun
        """

        # Run pipelines
        
        transformed_data = data.copy()
        
        transformed_data.columns = [col.replace(' ', '_') for col in transformed_data.columns]

        for pip in self.pipelines:

            self.log(log_message=f'Running pipeline: {pip.pipeline_name}', timestamp=True)

            transformed_data, analysis_log = pip.forward(transformed_data.copy())
            
            if len(analysis_log):
                
                self.log(log_message='ANALYSIS FLAGS:', timestamp=True)
            
            for l in analysis_log:
                
                self.log(log_message=f'  ** {l}', timestamp=False)

        return transformed_data

### ---------------------------------------- ###

def one_hot_encode_categoricals(data: pd.DataFrame, min_counts: int=3, min_repeated_vals: float=0.75) -> tuple[pd.DataFrame, list[str]]:
    
    """
    One-hot encodes categorical variables
    
    Parameters
    ----------
    data : pd.DataFrame
        Dataframe to process
    min_counts : int=3
        Minimum number of times a variable value should repeat to be considered a category
    min_repeated_vals : float=0.75
        Minimum amount of repeated values making up the variable for it to be considered categorical
    
    Returns
    -------
    onehot_encoded_data : pd.DataFrame
        Data with categorical variables one-hot encoded
    warnings : list[str]
        List of warning messages
    """
    
    # Reset index

    data = data.reset_index(drop=True)
    
    onehot_encoded_data = data.copy()
    
    # Parse columns
    
    warnings = []

    for col in data.columns:
        
        # Extract column data

        col_data = data[col]
        
        col_data = col_data.values
        
        if col_data.shape[0] == 0:
            
            continue

        # Get data type
        
        col_dtype = str(data[col].dtypes)
    
        # Check if data is categorical (need non float repeated values with count > min_counts make up >= min_repeated_vals of the data)
        
        if 'float' in col_dtype:
            
            categorical_check = False
        
        else:
        
            categories = {val : (col_data == val).sum() for val in set(col_data)}
            
            categorical_count = sum([n for c,n in categories.items() if n > min_counts])
            
            categorical_check = ((categorical_count / len(col_data)) >= min_repeated_vals)
        
        if categorical_check:
            
            try:
            
                onehot_encoded_col = pd.get_dummies(col_data)
                
                onehot_encoded_col.columns = (col + '_' + onehot_encoded_col.columns.astype(str)).tolist()
                
                onehot_encoded_data = pd.concat([onehot_encoded_data, onehot_encoded_col], axis=1)
                
                onehot_encoded_data = onehot_encoded_data.drop(col, axis=1)
            
            except:
                
                warnings.append(f'ONE HOT ENCODING FAILED for column {col}')
            
    return onehot_encoded_data, warnings

### ---------------------------------------- ###

def check_type_issues(data: pd.DataFrame, major_type_threshold: float=0.8) -> tuple[pd.DataFrame, list[str]]:
    
    """
    Checks if columns have mixed data types
    N.B. float and int are fine together
    
    Parameters
    ----------
    data : pd.DataFrame
        Dataframe to process
    major_type_threshold : float=0.8
        Minimum amount of values of a specific dtype to assign the dtype to the whole column
    
    Returns
    -------
    data : pd.DataFrame
        Original input data
    warnings : list[str]
        List of warning messages
    """

    # Reset index

    data = data.reset_index(drop=True)

    # Dict of type-checking functions

    types_to_check = {
        'int' : lambda val: (type(int(val)) == int) & (int(val) == float(val)),
        'float' : lambda val: type(float(val)) == float,
        'bool' : lambda val: val.lower() in ('true', 'false', 't', 'f', 'yes', 'no', 'y', 'n', 'on', 'off', '1', '0'),
        'datetime_fmt0' : lambda val: type(val) == datetime,
        'datetime_fmt1' : lambda val: type(datetime.strptime(val.replace('/', '-'), "%Y-%m-%d %H:%M:%S")) == datetime.datetime,
        'datetime_fmt2' : lambda val: type(datetime.strptime(val.replace('/', '-'), "%Y-%B-%d %H:%M:%S")) == datetime.datetime,
        'datetime_fmt3' : lambda val: type(datetime.strptime(val.replace('/', '-'), "%Y-%b-%d %H:%M:%S")) == datetime.datetime,
        'datetime_fmt4' : lambda val: type(datetime.strptime(val.replace('/', '-'), "%d-%m-%Y %H:%M:%S")) == datetime.datetime,
        'datetime_fmt5' : lambda val: type(datetime.strptime(val.replace('/', '-'), "%d-%B-%Y %H:%M:%S")) == datetime.datetime,
        'datetime_fmt6' : lambda val: type(datetime.strptime(val.replace('/', '-'), "%d-%b-%Y %H:%M:%S")) == datetime.datetime,
        'datetime_fmt7' : lambda val: type(datetime.strptime(val.replace('/', '-'), "%Y-%m-%d")) == datetime.datetime,
        'datetime_fmt8' : lambda val: type(datetime.strptime(val.replace('/', '-'), "%Y-%B-%d")) == datetime.datetime,
        'datetime_fmt9' : lambda val: type(datetime.strptime(val.replace('/', '-'), "%Y-%b-%d")) == datetime.datetime,
        'datetime_fmt10' : lambda val: type(datetime.strptime(val.replace('/', '-'), "%d-%m-%Y")) == datetime.datetime,
        'datetime_fmt11' : lambda val: type(datetime.strptime(val.replace('/', '-'), "%d-%B-%Y")) == datetime.datetime,
        'datetime_fmt12' : lambda val: type(datetime.strptime(val.replace('/', '-'), "%d-%b-%Y")) == datetime.datetime,
        'datetime_fmt13' : lambda val: type(datetime.strptime(val, "%H:%M:%S")) == datetime.datetime,
        'str' : lambda val: type(str(val)) == str
    }
    
    # Parse columns
    
    warnings = []
    
    data_dtypes = []

    for col in data.columns:
        
        # Extract column data

        col_data = data.loc[~ data[col].isna(), col]
        
        col_data_idx = col_data.index.values
        
        col_data = col_data.values
        
        if col_data.shape[0] == 0:
            
            data_dtypes.append(('unknown', False))
            
            continue

        # Check data types
        
        if str(data[col].dtypes) != 'object':
            
            col_dtype = str(data[col].dtypes)
        
        else:
        
            col_dtype = {type_name : 0 for type_name in types_to_check.keys()}
            
            d_types = []
            
            for d in col_data:
                
                assigned_dtype = 'unknown'
                
                for type_name,type_check in types_to_check.items():
                    
                    try:
                        
                        test = type_check(d)
                    
                    except:
                        
                        test = False
                        
                    if test:
                        
                        col_dtype[type_name] += 1
                        
                        assigned_dtype = type_name
                        
                        break
                    
                d_types.append(assigned_dtype)
    
            # Check if one data type is dominant
            # N.B. if both int and float are present, then int values are considered float
            
            present_types = [t for t,count in col_dtype.items() if count > 0]
            
            if 'float' in present_types and 'int' in present_types:
                
                col_dtype['float'] += col_dtype['int']
                
                col_dtype['int'] = 0
                
                present_types.remove('int')
            
            major_type = [t for t in present_types if (col_dtype[t] / len(col_data)) >= major_type_threshold]
            
            # Warn if more than one data type is present
            
            if len(present_types) == 1:
                
                col_dtype = present_types[0]
            
            elif len(present_types) > 1 and len(major_type):
                
                col_dtype = major_type[0]
                
                wrong_type_idx = [str(N + 1) for N,t in zip(col_data_idx, d_types) if t != col_dtype]
                
                warnings.append(f'WRONG VALUE TYPES in column "{col}" of type "{major_type}" at rows: [{", ".join(wrong_type_idx)}]')
            
            elif len(present_types) > 1:
                
                col_dtype = 'unknown'
                
                warnings.append(f'MIXED TYPES in column "{col}": [{", ".join(present_types)}]')
            
            else:
                
                col_dtype = 'unknown'
                
                warnings.append(f'UNKNOWN TYPE for column "{col}"')
            
        data_dtypes.append(col_dtype)
    
    # Add log of data_dtypes types

    dtypes_log = []
    
    dtypes_summary = pd.DataFrame({
        'column' : data.columns,
        'dtype' : data_dtypes,
    })
    
    for datatype in np.sort(np.unique(dtypes_summary['dtype'].values)):
        
        columns = dtypes_summary.loc[dtypes_summary['dtype'] == datatype, 'column'].values.astype(str)
        
        count = len(columns)
        
        dtypes_log.append(f'DTYPE "{datatype}" found in {count} columns: [{", ".join(columns)}]')
    
    warnings = dtypes_log + warnings
    
    return data, warnings

### ---------------------------------------- ###

def check_missing_values(data: pd.DataFrame, missing_threshold: float=0.1, max_missing: float=0.5, max_tests: int=25) -> tuple[pd.DataFrame, list[str]]:
    
    """
    Checks if columns have more than missing_threshold missing values
    If > missing_threshold and < max_missing, a warning is issued
    If > max_missing, data is discarded
    
    Parameters
    ----------
    data : pd.DataFrame
        Dataframe to process
    missing_threshold : float=0.1
        Minimum percentage of missing values to trigger a warning
    max_missing : float=0.5
        Minimum percentage of missing values to discard the column
    max_tests : int=25
        Maximum number of tests to perform to check for Missing Not At Random
    
    Returns
    -------
    cleaned_data : pd.DataFrame
        Data where columns with > max_missing missing values are discarded
    warnings : list[str]
        List of warning messages
    """
    
    # Reset index

    data = data.reset_index(drop=True)
    
    data_cleaned = data.copy()
    
    # Get possible test columns
    
    acceptable_test_dtypes = ['float', 'float32', 'float64', 'object', 'categorical', 'bool']
    
    test_cols = [col for col in data.columns if str(data[col].dtypes) in acceptable_test_dtypes]
    
    if len(test_cols) > max_tests:
        
        random_choice = np.random.default_rng(42)
        
        test_cols = random_choice.choice(test_cols, max_tests)
    
    # Parse columns
    
    warnings = []
    
    for col in data.columns:
        
        col_data = data[col]
        
        na_count = col_data.isna().sum() + (data[col].values == '').sum()
        
        col_data = col_data.values
        
        na_freq = na_count / data.shape[0]
        
        if na_freq >= max_missing:
            
            data_cleaned = data_cleaned.drop(col, axis=1)
            
            warnings.append(f'MISSING VALUES in column "{col}" make up {(100 * na_freq):.3f}% of all values. Column was dropped')
        
        elif na_freq >= missing_threshold:
            
            # Check why data is missing not at random (MNAR)
            
            actual_tests, mnar_proof = 0, 0
            
            p_thr = 0.05
            
            for tc in test_cols:
                
                # Testing if the data in test_column has the same distribution when using all values or the ones masked as not missing in col 
                
                test_data_1 = data.loc[:, tc].dropna()
                
                test_data_2 = data.loc[~ data[col].isna(), tc].dropna()
                
                if len(test_data_1) < 3 or len(test_data_2) < 3:
                    
                    continue
                
                if 'float' in str(test_data_1.dtypes):
                    
                    # Continuous data
                    
                    try:

                        pval = mannwhitneyu(test_data_1, test_data_2, alternative='two-sided').pvalue

                        actual_tests += 1

                    except:

                        pval = 1
                
                else:
                    
                    # Categorical data
                    
                    unique_categories = np.unique(test_data_1)

                    test_data_1, test_data_2 = test_data_1.to_list(), test_data_2.to_list()
                    
                    test_data_1_freq = [test_data_1.count(c) for c in unique_categories]
                    
                    test_data_2_freq = [test_data_2.count(c) for c in unique_categories]
                    
                    test_data_2_freq = [round(c * sum(test_data_1_freq) / sum(test_data_2_freq)) for c in test_data_2_freq]
                    
                    try:
                    
                        pval = chisquare(test_data_2_freq, test_data_1_freq).pvalue

                        actual_tests += 1
            
                    except:
                        
                        pval = 1.
                
                if pval < p_thr:
                    
                    mnar_proof += 1

            mnar_p = (actual_tests - mnar_proof) / actual_tests
            
            if mnar_p < 0.05:
                
                warnings.append(f'MISSING VALUES in column "{col}" make up {(100 * na_freq):.3f}% of all values. Suspected MNAR')
            
            else:
            
                warnings.append(f'MISSING VALUES in column "{col}" make up {(100 * na_freq):.3f}% of all values.')
        
        else:
            
            continue
    
    return data_cleaned, warnings

### ---------------------------------------- ###

def outlier_detection(data: pd.DataFrame, remove_outliers: bool=False, mahalanobis_only: bool=False, multihit_required: bool=True) -> tuple[pd.DataFrame, list[str]]:
    
    """
    Checks for outliers using several methods
    
    Parameters
    ----------
    data : pd.DataFrame
        Dataframe to process
    remove_outliers : bool=False
        Set to True to remove samples with outlier values
    mahalanobis_only : bool=True
        Set to True to remove outliers found using Mahalanobis distance
        Use in with multihit_required=True to further subset
    multihit_required : bool=True
        Set to True to require more than one metric to be used to mark a sample as outlier before removal
    
    Returns
    -------
    data_cleaned : pd.DataFrame
        Original input data or data with outliers removed
    warnings : list[str]
        List of warning messages
    """
    
    # Reset index

    data = data.reset_index(drop=True)
    
    data_cleaned = data.copy()
    
    # Get float columns
    
    acceptable_test_dtypes = ['float', 'float32', 'float64']
    
    test_cols = [col for col in data.columns if str(data[col].dtypes) in acceptable_test_dtypes]
    
    # Make sure none of the test_cols is actually categorical
    
    cols_to_remove = []
    
    for col in test_cols:
        
        col_data = data.loc[~ data[col].isna(), col].values
    
        categories = {val : (col_data == val).sum() for val in set(col_data)}
    
        categorical_count = sum([n for c,n in categories.items() if n > 3])
    
        categorical_check = ((categorical_count / len(col_data)) >= 0.75)
        
        if categorical_check:
            
            cols_to_remove.append(col)
    
    test_cols = [col for col in test_cols if col not in cols_to_remove]
    
    # Check outliers in individual columns
    
    warnings = []
    
    outlier_rows = []
    
    for col in test_cols:
        
        # Extract column data

        col_data = data.loc[~ data[col].isna(), col]
        
        col_data_idx = col_data.index.values
        
        col_data = col_data.values
        
        if col_data.shape[0] == 0:
            
            continue
        
        # IQR
        
        q2 = np.median(col_data)
        
        q1 = np.median(col_data[col_data < q2])
        
        q3 = np.median(col_data[col_data > q2])
        
        iqr = q3 - q1
        
        thr_low = q1 - 1.5 * iqr
        
        thr_high = q3 + 1.5 * iqr
        
        outliers = col_data_idx[np.where((col_data < thr_low) | (col_data > thr_high))[0]]
        
        if len(outliers):
            
            outlier_rows += outliers.tolist()
            
            warnings.append(f'OUTLIERS FOUND using IQR for columns "{col}" at rows: [{", ".join(outliers.astype(str))}]')
        
        # Zscore
        
        zscaled = (col_data.copy() - col_data.mean()) / col_data.std()
        
        z_thr = 3
        
        outliers = col_data_idx[np.where((zscaled < - z_thr) | (zscaled > z_thr))[0]]
        
        if len(outliers):
            
            outlier_rows += outliers.tolist()
        
            warnings.append(f'OUTLIERS FOUND using Zscore for columns "{col}" at rows: [{", ".join(outliers.astype(str))}]')
        
        # MAD Zscore
        
        col_median = np.median(col_data)
        
        mad = np.median(np.abs(col_data - col_median))
        
        mad_zscaled = 0.6745 * (col_data.copy() - col_median) / mad
        
        z_thr = 3
        
        outliers = col_data_idx[np.where((mad_zscaled < - z_thr) | (mad_zscaled > z_thr))[0]]
        
        if len(outliers):
            
            outlier_rows += outliers.tolist()
        
            warnings.append(f'OUTLIERS FOUND using MAD Zscore for columns "{col}" at rows: [{", ".join(outliers.astype(str))}]')
        
    # Check outlier using all columns with Mahalanobis distance
    # N.B.  The squared Mahalanobis distance follows a chi-squared distribution
    
    data_sub = data[test_cols]
    
    data_sub = data_sub.loc[:, data_sub.isna().sum(axis=0) == 0] # Only use columns without missing values

    data_sub_mean = data_sub.mean(axis=0)
    
    data_sub_cov_matrix = np.cov(data_sub, rowvar=False)
    
    data_sub_inverse_cov_matrix = np.linalg.inv(data_sub_cov_matrix)
    
    m_distances = [distance.mahalanobis(row.values, data_sub_mean, data_sub_inverse_cov_matrix) for _,row in data_sub.iterrows()]
    
    m_distances_squared = np.power(m_distances, 2)
    
    m_thr = chi2.ppf(0.95, len(test_cols))
    
    outliers = col_data_idx[np.where(m_distances_squared > m_thr)[0]]

    mahalanobis_outliers = []
    
    if len(outliers):
        
        outlier_rows += outliers.tolist()

        mahalanobis_outliers += outliers.tolist()
    
        warnings.append(f'OUTLIERS FOUND using Mahalanobis distance at rows: [{", ".join(outliers.astype(str))}]')
    
    # Remove outliers
    
    if remove_outliers:

        if mahalanobis_outliers and multihit_required:

            outlier_rows = [orow for orow in outlier_rows if outlier_rows.count(orow) > 1 and orow in mahalanobis_outliers]
        
        elif mahalanobis_outliers and not multihit_required:
            
            outlier_rows = mahalanobis_outliers

        elif not mahalanobis_outliers and multihit_required:

            outlier_rows = [orow for orow in outlier_rows if outlier_rows.count(orow) > 1]

        else:

            pass
        
        outlier_rows = np.sort(np.unique(outlier_rows))

        data_cleaned = data.loc[~ data.index.isin(outlier_rows),]

        warnings.append(f'OUTLIERS ROWS REMOVED: {len(outlier_rows)} [{", ".join(outlier_rows.astype(str))}]')
    
    return data_cleaned, warnings

### ---------------------------------------- ###

def find_correlated_variables(data: pd.DataFrame, corr_thr: float=0.75) -> tuple[pd.DataFrame, list[str]]:
    
    """
    Check for correlated variables
    
    Parameters
    ----------
    data : pd.DataFrame
        Dataframe to process
    corr_thr : float=0.75
        Minimum coefficient of determination value to mark variables as correlated
    
    Returns
    -------
    data : pd.DataFrame
        Original input data
    warnings : list[str]
        List of warning messages
    """
    
    # Reset index

    data = data.reset_index(drop=True)
    
    # Get acceptable test columns
    
    acceptable_test_dtypes = ['float', 'float32', 'float64', 'categorical', 'bool']
    
    test_cols = [col for col in data.columns if str(data[col].dtypes) in acceptable_test_dtypes]
    
    # Correlation matrices
    
    warnings = []
    
    for corr_method in ['pearson', 'spearman', 'kendall']:
        
        # Correlation matrix
        
        corr_matrix = data[test_cols].corr(method=corr_method, numeric_only=False)
        
        corr_matrix_squared = corr_matrix ** 2
        
        # Cleanup pairs
        
        a, b = np.where(corr_matrix_squared >= corr_thr)
        
        a, b = a[a != b], b[b != a] # Removing identities
        
        a, b = a[a < b], b[b > a] # Removing symmetric
        
        # Cluster
        
        pairs = np.array([a, b]).T
        
        seeds = np.unique(pairs[~ np.isin(pairs[:, 0], pairs[:, 1]), 0])
        
        clusters, clustered = [], set()
        
        for s in seeds:
            
            if s in clustered:
                
                continue
            
            else:
                
                clustered.update({s})
            
            new_cluster = {s}
            
            toggle = True
            
            while toggle:
                
                elements_rows = np.where(np.isin(pairs, list(new_cluster)))[0]
                
                new_cluster_updated = set(pairs[elements_rows,].ravel())
                
                new_cluster_updated.update(new_cluster)
                
                if new_cluster_updated == new_cluster:
                    
                    toggle = False
                
                new_cluster = new_cluster_updated
                
            clusters.append(new_cluster)
            
            clustered.update(new_cluster)
        
        for cl in clusters:
            
            cl_elements = [f'"{test_cols[c_idx]}"' for c_idx in cl]
            
            warnings.append(f'CORRELATED VARIABLES FOUND using {corr_method} correlation: [{", ".join(cl_elements)}]')
    
    return data, warnings

### ---------------------------------------- ###

def check_distribution(data: pd.DataFrame, apply_transform: bool=False) -> tuple[pd.DataFrame, list[str]]:
    
    """
    Check data for normality and attempts to correct data
    
    Parameters
    ----------
    data : pd.DataFrame
        Dataframe to process
    apply_transform : bool=False
        Set to True to apply best proposed transformations for improving normality
    
    Returns
    -------
    data : pd.DataFrame
        Original input data or transformed data
    warnings : list[str]
        List of warning messages
    """
    
    # Reset index

    data = data.reset_index(drop=True)
    
    transformed_data = data.copy()
    
    # Get float columns
    
    acceptable_test_dtypes = ['float', 'float32', 'float64']
    
    test_cols = [col for col in data.columns if str(data[col].dtypes) in acceptable_test_dtypes]
    
    # Make sure none of the test_cols is actually categorical
    
    cols_to_remove = []
    
    for col in test_cols:
        
        col_data = data.loc[~ data[col].isna(), col].values
    
        categories = {val : (col_data == val).sum() for val in set(col_data)}
    
        categorical_count = sum([n for c,n in categories.items() if n > 3])
    
        categorical_check = ((categorical_count / len(col_data)) >= 0.75)
        
        if categorical_check:
            
            cols_to_remove.append(col)
    
    test_cols = [col for col in test_cols if col not in cols_to_remove]
    
    # Parse columns
    
    warnings = []
    
    transformations = {
        'none' : lambda dt: dt,
        'log1p' : lambda dt: np.log1p(dt),
        'sqrt' : lambda dt: np.sqrt(dt),
        'cbrt' : lambda dt: np.cbrt(dt),
        'boxcox' : lambda dt: boxcox(dt)[0]
    }
    
    params_fail_check = {
        'skweness' : lambda val: abs(val) > 0.5,
        'kurtosis' : lambda val: abs(val) > 1,
        'normality_p' : lambda val: val < 0.05
    }
    
    for col in test_cols:
        
        # Extract column data

        col_data = data.loc[~ data[col].isna(), col]
        
        if col_data.shape[0] == 0:
            
            continue
    
        # Test transformations
        
        distribution_parameters = {param : {} for param in ['skweness', 'kurtosis', 'normality_p']}
        
        with catch_warnings():
            
            simplefilter("ignore", category=RuntimeWarning)
        
            for tr_name,tr_fun in transformations.items():
                
                try:
                    
                    # Transform
                    
                    col_data_tr = tr_fun(col_data)
                    
                    # Check skewness and kurtosis
                    
                    distribution_parameters['skweness'][tr_name] = col_data_tr.skew()
                    
                    distribution_parameters['kurtosis'][tr_name] = col_data_tr.kurtosis()
                    
                    distribution_parameters['normality_p'][tr_name] = normaltest(col_data_tr, axis=0, nan_policy='omit').pvalue
                
                except:
                    
                    pass
        
        # Suggest transformations
        
        for param, values in distribution_parameters.items():
            
            # Check if within thresholds
            
            values_fail = {tr_name : params_fail_check[param](tr_value) for tr_name,tr_value in values.items()}
            
            helpful_tr = [tr_name for tr_name,vf in values_fail.items() if not vf and tr_name != 'none']
            
            if values_fail['none']:
                
                # Raw value is not great
                
                if len(helpful_tr):
                    
                    # Some transformations can help
                
                    warnings.append(f'DISTRIBUTION ISSUE for column "{col}": {param}. The following transformations can help: [{", ".join(helpful_tr)}]')
                
                else:
                    
                    # Transformations cannot help
                    
                    warnings.append(f'DISTRIBUTION ISSUE for column "{col}": {param}')
            
            if apply_transform and param == 'normality_p' and values_fail['none']:
                
                with catch_warnings():
                    
                    simplefilter("ignore", category=RuntimeWarning)
                    
                    best_tr = max(values, key=values.get)

                    if best_tr in helpful_tr:
                
                        col_data_tr = transformations[best_tr](col_data)
                        
                        transformed_data.loc[:, col] = col_data_tr
                        
                        warnings.append(f'APPLIED TRANSFORMATION "{best_tr}" for column "{col}". Normality pvalue improved from {distribution_parameters["normality_p"]["none"]:.2e} to {distribution_parameters["normality_p"][best_tr]:.2e}')
    
    return transformed_data, warnings

### ---------------------------------------- ###

if __name__ == '__main__':

    # Additional imports

    from sys import argv
    from sys import exit as sys_exit
    from os.path import exists
    
    # Parse cli args

    data_path, sheet, one_hot, remove_outliers, apply_transform = parse_args()

    # Run data QC

    if exists(data_path):
        
        # Load data
        
        if any(data_path.endswith(suffix) for suffix in ['.tsv', '.csv', '.txt']):
        
            sep = '\t' if data_path.endswith('.tsv') else ',' if data_path.endswith('.csv') else ','
            
            data = pd.read_csv(data_path, sep='\t')
        
        elif data_path.endswith('.xlsx'):
            
            try:
                
                sheet = int(sheet)
            
            except:
                
                pass
            
            data = pd.read_excel(data_path, sheet_name=sheet)
        
        else:
            
            print('ERROR: unrecognized file extension')
            
            sys_exit()

        # Init pipeline
        
        manager = pipelines_manager()
        
        # Add pipelines
        
        if one_hot:
        
            manager.add_pipeline(
                pipeline(
                    pipeline_name='One-hot encoding',
                    pipeline_description='Finds categorical columns and one-hot encodes them',
                    analysis_fun=lambda dt: one_hot_encode_categoricals(data=dt, min_counts=3, min_repeated_vals=0.75)
                )
            )
        
        manager.add_pipeline(
            pipeline(
                pipeline_name='Check dtypes',
                pipeline_description='Checks and reports issues with columns dtypes',
                analysis_fun=lambda dt: check_type_issues(data=dt, major_type_threshold=0.8)
            )
        )
        
        manager.add_pipeline(
            pipeline(
                pipeline_name='Check missing values',
                pipeline_description='Checks for missing values, removes columns with too many, and checks if missingness is not random',
                analysis_fun=lambda dt: check_missing_values(data=dt, missing_threshold=0.1, max_missing=0.5, max_tests=25)
            )
        )
        
        manager.add_pipeline(
            pipeline(
                pipeline_name='Outlier detection',
                pipeline_description='Detects outliers using several methods and removes them if desired',
                analysis_fun=lambda dt: outlier_detection(data=dt, remove_outliers=remove_outliers, mahalanobis_only=True, multihit_required=True)
            )
        )
        
        manager.add_pipeline(
            pipeline(
                pipeline_name='Variables correlation',
                pipeline_description='Finds correlated variables using several methods',
                analysis_fun=lambda dt: find_correlated_variables(data=dt, corr_thr=0.75)
            )
        )
        
        manager.add_pipeline(
            pipeline(
                pipeline_name='Check distribution',
                pipeline_description='Checks skewness, kurtosis, and normality and suggests (or applies) transformations to resolve issues',
                analysis_fun=lambda dt: check_distribution(data=dt, apply_transform=apply_transform)
            )
        )
        
        # Run exploration
        
        transformed_data = manager.explore(data)
        
        transformed_data.to_csv('transformed_data.tsv', sep='\t', index=False)
