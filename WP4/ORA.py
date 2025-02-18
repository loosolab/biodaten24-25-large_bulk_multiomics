#!/usr/bin/env python3
import scipy.stats as stats
import pandas as pd

class ORA:

    def __init__(self):
        pass

    def categorical(self, meta_column, cl_samples):
        # N = Total population (e.g., total samples)
        # M = Total number of successes (e.g., total males)
        # n = Sample size (e.g., size of the cluster)
        # k = Number of successes in the sample (e.g., number of males in the cluster)

        N = len(meta_column)
        n = len(cl_samples)
        pvals_dict = {}
        meta_dict = meta_column.to_dict()
        category_counts = meta_column.value_counts().to_dict()
        cluster_categories = [meta_dict[sample] for sample in cl_samples if sample in meta_dict]
        cluster_counts = pd.Series(cluster_categories).value_counts().to_dict()
        for category, M in category_counts.items():
            k = cluster_counts.get(category, 0)
            p_value = stats.hypergeom.sf(k - 1, N, M, n)
            pvals_dict[category] = p_value
        return pvals_dict