import random
import numpy as np
import pandas as pd

class BTSing:
    def __init__(self):
        pass

    def categorical(self, meta_column, cl_samples, n_iter=1000):
        metacolumn_dict = meta_column.to_dict()
        unique_vars = meta_column.unique()
        temp_counts = {var: sum(1 for s in cl_samples if metacolumn_dict.get(s) == var) for var in unique_vars}
        perm_counts = {var: [] for var in unique_vars}
        for _ in range(n_iter):
            sample_size = min(len(cl_samples), len(meta_column))
            sampled = np.random.choice(meta_column.index, sample_size, replace=False)
            for var in unique_vars:
                perm_counts[var].append(sum(1 for s in sampled if metacolumn_dict.get(s) == var))
        pvals_dict = {}
        for var in unique_vars:
            perm_counts[var].sort(reverse=True)
            pvals_dict[var] = sum(1 for count in perm_counts[var] if count >= temp_counts[var]) / n_iter

        return pvals_dict
        
    def numerical(self, meta_column, cl_samples, n_iter=1000):
        cluster_values = meta_column.loc[cl_samples].dropna().values
        non_cluster_values = meta_column.drop(cl_samples, errors='ignore').dropna().values
        if len(cluster_values) == 0 or len(non_cluster_values) == 0:
            print ("Cluster or non-cluster samples are empty. Check input data.")
            pval = "ValueError"
            return pval
        observed_diff = abs(np.mean(cluster_values) - np.mean(non_cluster_values))
        combined_values = np.concatenate([cluster_values, non_cluster_values])
        boot_diffs = []
        for _ in range(n_iter):
            np.random.shuffle(combined_values)
            random_cluster = np.random.choice(combined_values, size=len(cluster_values), replace=False)
            random_non_cluster = np.random.choice(combined_values, size=len(non_cluster_values), replace=False)
            boot_diffs.append(abs(np.mean(random_cluster) - np.mean(random_non_cluster)))
        boot_diffs = np.array(boot_diffs)
        pval = np.sum(boot_diffs >= observed_diff) / n_iter
        return pval
    
    
