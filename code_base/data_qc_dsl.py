#!/usr/bin/env python3

"""
Domain-Specific Language for general data quality control
"""

### IMPORTS -------------------------------- ###

import pandas as pd

from collections.abc import Callable
from datetime import datetime

### CLASSES AND FUNCTIONS ------------------ ###

def help_fun():
    
    print("""
    Domain-Specific Language for general data quality control
    
    USAGE
    -----
        python data_qc_dsl.py [args]
    
    ARGS
    ----
        --data_path
            Path to data file
            Required
        --sheet
            Sheet of xlsx file to be used
            Omit if file is text
""")

### ---------------------------------------- ###

def parse_args() -> tuple[str]:
    
    """
    Parses CLI arguments
    
    Returns
    -------
    data_path : str
        Path to data file
    sheet : str
        Sheet of xlsx file to be used
    """
    
    # Path to data (tsv, csv, txt, or xlsx table)
    
    data_path = argv[argv.index('--data_path') + 1] if '--data_path' in argv else ''
    
    # Sheet for xlsx file
    
    sheet = argv[argv.index('--sheet') + 1] if '--sheet' in argv else '0'
    
    return data_path, sheet

### ---------------------------------------- ###

class qc_rule:

    """
    Main data QC-check class
    
    Parameters
    ----------
    rule_name : str
        Name used by the rule
    qc_fun : Callable
        Function to apply to the dataset
    error_message : str
        Description of errors identified by the rule
    
    Methods
    -------
    forward(data: pd.DataFrame)
        Applies the analysis_fun to data
    """

    def __init__(self, rule_name: str, qc_fun: Callable, error_message: str) -> None:
        
        """
        Class init
        
        Parameters
        ----------
        rule_name : str
            Name used by the rule
        qc_fun : Callable
            Function to apply to the dataset
        error_message : str
            Description of errors identified by the rule
        """

        self.rule_name = rule_name

        self.qc_fun = qc_fun

        self.error_message = error_message

    ### ------------------------------------ ###

    def forward(self, data: pd.DataFrame) -> dict[str]:
        
        """
        Applies the qc_fun to data
        
        Parameters
        ----------
        data : pd.DataFrame
            Dataframe to process
        
        Returns
        -------
        report : dict[str]
            Dictionary containing a report of the QC check
        """

        # Find issues
        
        issues = self.qc_fun(data)
        
        # Generate report
        
        report = {
            'rule_name' : self.rule_name,
            'issue_type' : self.error_message,
            'issues' : len(issues),
            'details' : issues
        }

        return report

### ---------------------------------------- ###

class quality_inspector:

    """
    Class handling QC rules execution
    
    Methods
    -------
    log(log_message: str, timestamp: bool=True)
        Prints a log message to stdout
    add_rule(new_rule: Callable)
        Adds a qc_rule class object
    inspect(data: pd.DataFrame)
        Runs QC rules sequentially on the data
    """

    def __init__(self) -> None:
        
        """
        Class init
        """

        self.qc_rule_set = []
        
        self.log(log_message='Quality inspector initialized', timestamp=True)

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

    def add_rule(self, new_rule: Callable) -> None:
        
        """
        Adds a qc_rule class object
        
        Parameters
        ----------
        new_rule : Callable
            New qc_rule class object to add
        """

        self.qc_rule_set.append(new_rule)

        self.log(log_message=f'Added new rule: {new_rule.rule_name}', timestamp=True)

    ### ------------------------------------ ###

    def inspect(self, data: pd.DataFrame) -> list[dict]:
        
        """
        Runs QC rules sequentially on the data
        
        Parameters
        ----------
        data : pd.DataFrame
            Dataframe to process
        
        Returns
        -------
        reports : list[dict]
            List of QC reports
        """
        
        # Inspect
        
        self.log(log_message='Inspecting data', timestamp=True)
        
        reports = []

        for rule in self.qc_rule_set:

            reports.append(rule.forward(data))
        
        # Format report
        
        self.log(log_message='Report ready:', timestamp=True)
        
        for r in reports:
            
            self.log(log_message='-' * 40, timestamp=False)
            self.log(log_message=f'### RULE: {r["rule_name"]}', timestamp=False)
            self.log(log_message=f'  * ISSUE: {r["issue_type"]}', timestamp=False)
            self.log(log_message=f'  * ISSUE COUNT: {r["issues"]}', timestamp=False)
            
            for detail in r['details']:
                
                self.log(log_message=f'    \u2717 {detail}', timestamp=False)
        
            self.log(log_message='-' * 40, timestamp=False)
    
        return reports

### ---------------------------------------- ###

def check_column_data_type_consistency(data: pd.DataFrame, major_type_threshold: float=0.8) -> list[str]:
    
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

    for col in data.columns:
        
        # Extract column data

        col_data = data.loc[~ data[col].isna(), col]
        
        col_data_idx = col_data.index.values
        
        col_data = col_data.values

        # Check data types
        
        col_types = {type_name : 0 for type_name in types_to_check.keys()}
        
        d_types = []
        
        for d in col_data:
            
            assigned_dtype = 'unknown'
            
            for type_name,type_check in types_to_check.items():
                
                try:
                    
                    test = type_check(d)
                
                except:
                    
                    test = False
                    
                if test:
                    
                    col_types[type_name] += 1
                    
                    assigned_dtype = type_name
                    
                    break
                
            d_types.append(assigned_dtype)

        # Check if one data type is dominant
        # N.B. if both int and float are present, then int values are considered float
        
        present_types = [t for t,count in col_types.items() if count > 0]
        
        if 'float' in present_types and 'int' in present_types:
            
            col_types['float'] += col_types['int']
            
            col_types['int'] = 0
            
            present_types.remove('int')
        
        major_type = [t for t in present_types if (col_types[t] / len(col_data)) >= major_type_threshold]
        
        # Warn if more than one data type is present
        
        if len(present_types) > 1 and len(major_type):
            
            major_type = major_type[0]
            
            wrong_type_idx = [str(N + 1) for N,t in zip(col_data_idx, d_types) if t != major_type]
            
            warnings.append(f'WRONG VALUE TYPES in column "{col}" of type "{major_type}" at rows: [{", ".join(wrong_type_idx)}]')
        
        elif len(present_types) > 1:
            
            warnings.append(f'MIXED TYPES in column "{col}": [{", ".join(present_types)}]')
        
        else:
            
            continue
        
    return warnings

### ---------------------------------------- ###

def check_missing_header_values(data: pd.DataFrame) -> list[str]:
    
    """
    Checks if columns header values are missing
    
    Parameters
    ----------
    data : pd.DataFrame
        Dataframe to process
    
    Returns
    -------
    warnings : list[str]
        List of warning messages
    """
    
    # Parse columns
    
    warnings = []
    
    for N,col in enumerate(data.columns):
        
        if col.startswith('Unnamed:') or col == '':
            
            warnings.append(f'MISSING HEADER VALUE at column {N + 1}')
        
    return warnings

### ---------------------------------------- ###

def check_missing_values(data: pd.DataFrame, missing_threshold: float=0.1):
    
    """
    Checks if columns have more than missing_threshold missing values
    
    Parameters
    ----------
    data : pd.DataFrame
        Dataframe to process
    missing_threshold : float=0.1
        Minimum percentage of missing values to trigger a warning
    
    Returns
    -------
    warnings : list[str]
        List of warning messages
    """
    
    # Parse columns
    
    warnings = []
    
    for col in data.columns:
        
        na_count = data[col].isna().sum() + (data[col].values == '').sum()
        
        na_freq = na_count / data.shape[0]
        
        if na_freq >= missing_threshold:
            
            warnings.append(f'MISSING VALUES in column "{col}" make up {(100 * na_freq):.3f}% of all values')
    
    return warnings

### ---------------------------------------- ###

def check_rare_categorical_values(data: pd.DataFrame, min_counts: int=3, min_repeated_vals: float=0.75) -> list[str]:
    
    """
    Checks if categorical columns have rare categories and check if they're likely to be due to typos
    
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
    warnings : list[str]
        List of warning messages
    """
    
    # Parse columns
    
    warnings = []
    
    for col in data.columns:
        
        # Extract column data

        col_data = data.loc[~ data[col].isna(), col].values
        
        if col_data.shape[0] == 0:
            
            continue
        
        # Check if data is categorical
        
        categories = {val : (col_data == val).sum() for val in set(col_data)}
        
        # Skip if repeated values with count > min_counts make up less than min_repeated_vals of the data
        
        categorical_count = sum([n for c,n in categories.items() if n > min_counts])
        
        if (categorical_count / len(col_data)) < min_repeated_vals:
            
            continue
        
        # Skip if there's too many categories proportionally to the number of data points
        
        max_categories = int(0.1 * len(col_data))
        
        if len(categories) > max_categories:
            
            continue
        
        # Define a threshold for rarity using kneedle method
        
        rarity_threshold = []
        
        counts = list(categories.values())
        counts.sort(reverse=True)
        
        x0, x1 = 0, len(counts)
        y0, y1 = max(counts), min(counts)
        gradient = (y1 - y0) / (x1 - x0)
        intercept = y0
        difference_vector = [(gradient * x + intercept) - y for x,y in enumerate(counts)]
        cutoff_index = difference_vector.index(max(difference_vector))
        
        #min_frequency_drop = 0.25
        #drop_check = ((counts[cutoff_index - 1] - counts[cutoff_index]) / counts[cutoff_index - 1]) > min_frequency_drop
        #rarity_threshold = counts[cutoff_index] if drop_check else counts[-1] - 1
        rarity_threshold = counts[cutoff_index]
        
        # Flag rare categories
        
        rare_categories = [c for c,n in categories.items() if n < rarity_threshold]
        
        for rc in rare_categories:
            
            # Check for categories with similar names if rc is a string
            # i.e. check if category is a typo
            
            if type(rc) == str:
            
                similarly_named_categories = [(common, sum([a == b for a,b in zip(rc, common)]) / min(len(rc), len(common))) for common in categories.keys() if common not in rare_categories and type(common) == str]
                
                similarly_named_categories = [(category, score) for category,score in similarly_named_categories if score >= 0.5]
            
            else:
                
                similarly_named_categories = []
            
            if len(similarly_named_categories):
                
                best_matching_value = max(similarly_named_categories, key=lambda c: c[1])[0]
                
                warnings.append(f'RARE CATEGORICAL in column "{col}" similar to value "{best_matching_value}": "{rc}"')
            
            else:
                
                warnings.append(f'RARE CATEGORICAL in column "{col}": "{rc}"')
    
    return warnings

### ---------------------------------------- ###

if __name__ == '__main__':

    # User provided file path

    from sys import argv
    from sys import exit as sys_exit
    from os.path import exists

    # Parse cli args

    data_path, sheet = parse_args()

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

        # Init inspector
        
        inspector = quality_inspector()
        
        # Add tools
        
        inspector.add_rule(
            qc_rule(
                'Inconsistent data types detection',
                lambda dt: check_column_data_type_consistency(data=dt, major_type_threshold=0.8),
                'WARNING: inconsistent data types detected'
            )
        )
        
        inspector.add_rule(
            qc_rule(
                'Missing header names',
                check_missing_header_values,
                'WARNING: header incomplete'
            )
        )
        
        inspector.add_rule(
            qc_rule(
                'Missing values',
                lambda dt: check_missing_values(data=dt, missing_threshold=0.1),
                'WARNING: columns have more missing values than tolerated'
            )
        )
        
        inspector.add_rule(
            qc_rule(
                'Rare categorical values',
                lambda dt: check_rare_categorical_values(data=dt, min_counts=3, min_repeated_vals=0.75),
                'WARNING: some columns may have categorical values with typos'
            )
        )
        
        # Check data
        
        _ = inspector.inspect(data)
