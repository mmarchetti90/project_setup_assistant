#!/usr/bin/env python3

### ---------------------------------------- ###

class fit_glm():
    
    """
    Class for fitting a Generalized Linear Model to a set of variables.

    Parameters
    ----------
    endog : pd.DataFrame
        Data frame of dipendent variables with shape NxM, with N being data points and M individual
        variables.
        The index of the data frame must match those found in exog.
    exog : pd.DataFrame
        Data frame of indipendent variables with shape NxM, with N being data points and M
        individual variables.
        The index of the data frame must match those found in endog.

    Attributes
    ----------
    endog : pd.DataFrame
        Data frame of dipendent variables.
    exog : pd.DataFrame
        Data frame of indipendent variables.

    Methods
    -------
    fit()
        Fits a GLM to each endog variable.
    """
    
    def __init__(self, endog, exog):
        
        endog.index = [idx.replace(' ', '_') for idx in endog.index.values]
        endog.columns = [col.replace(' ', '_') for col in endog.columns.values]
        self.endog = endog.copy()
        
        exog.index = [idx.replace(' ', '_') for idx in exog.index.values]
        self.exog = exog.copy()
        
    ### ------------------------------------ ###
    
    def fit(self, formula_right, model_family=None, save_results=True, output_prefix='glm_results'):
        
        """
        Fits a GLM to each endog variable.
        Returns a summary results data frame.

        Parameters
        ----------
        formula_right : str
            Right part of the formula to be used during fitting.
        model_family : function
            Family function from statsmodels.genmod.families.family to be passed to
            statsmodels.genmod.generalized_linear_model.GLM
            Default=None (i.e. statsmodels.genmod.families.family.Gaussian)
        save_results  bool
            If True, saves results and fitting log as text files.
            Default=True
        output_prefix : str
            Prefix for output files.
            Default='glm_results'
        """
        
        # Init results table and log
        
        results = []
        results_header = ['variable', 'parameter', 'value', 'coef', 'pval']
        log = ['#' * 40,
               'Fitting models with the following formula:',
               f'~ {formula_right}',
               '#' * 40]
        
        # Fit models for each variable
        
        for var in self.endog.columns.values:
            
            endog_sub = self.endog[var].copy()
            
            if endog_sub.sum() == 0:
                
                continue
            
            model_dat = pd.merge(endog_sub, self.exog.copy(), how='inner', left_index=True, right_index=True)
            
            formula = f'{var} ~ {formula_right}'
            
            model = GLM.from_formula(formula, model_dat, family=model_family)
            
            res = model.fit()
            
            # Store results
            
            coefs, pvals = res.params, res.pvalues
            coefs.name, pvals.name = 'coef', 'pval'
            
            for idx,(coef,pval) in pd.merge(coefs, pvals, how='inner', left_index=True, right_index=True).iterrows():
                
                if idx == 'Intercept':
                    
                    param, value = 'Intercept', ''
                    
                else:
                    
                    idx = idx.split(':')
                    
                    param = [i if '[' not in i else
                             i[:i.index('[')].replace('C(', '').replace(')', '')
                             for i in idx]
                    
                    value = [i[i.index('['):].replace('[T.', '').replace(']', '')
                             for i in idx
                             if '[' in i]
                
                    param, value = '_*_'.join(param), '_*_'.join(value)
                
                results.append([var, param, value, coef, pval])
            
            # Update log
            
            log.append('#' * 40)
            log.append(f'\nFitting effects on {var}\n')
            log.append(str(res.summary()))
            log.append('')
            
        results = pd.DataFrame(results, columns=results_header)
        
        ### Add baseline info for categorical variables

        results['baseline'] = np.repeat('', results.shape[0])

        categorical_vars = np.unique(results.loc[results['value'] != '', 'parameter'])

        for cat_var in categorical_vars:
            
            baseline = []
            
            for cv in cat_var.split('_*_'):
            
                b = np.sort(np.unique(self.exog[cv].values))[0]
                
                baseline.append(str(b))
            
            baseline = '_*_'.join(baseline)
            
            results.loc[results['parameter'] == cat_var, 'baseline'] = baseline
        
        results = results[['variable', 'parameter', 'baseline', 'value', 'coef', 'pval']]
        
        ### Adjust pvalues

        results['padj'] = fdrcorrection(results['pval'].values, alpha=0.05, is_sorted=False)[1]
            
        results.sort_values(by='pval', ascending=True, inplace=True)

        ### Save data

        if save_results:
            
            results.to_csv(f'{output_prefix}.tsv', sep='\t', index=False, header=True)
    
            with open(f'{output_prefix}.log', 'w') as log_out:
                    
                log_out.write('\n'.join(log))
        
        return results

### ------------------MAIN------------------ ###

import numpy as np
import pandas as pd

from statsmodels.genmod.generalized_linear_model import GLM
from statsmodels.stats.multitest import fdrcorrection
