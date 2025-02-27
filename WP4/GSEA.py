#!/usr/bin/env python3
import os
import pandas as pd
import gseapy as gp
from pathlib import Path

class GSEA:
    #path management
    parent_dir = Path(__file__).parent
    os.makedirs(f"{parent_dir}/gsea/input", exist_ok=True)
    os.makedirs(f"{parent_dir}/gsea/output", exist_ok=True)
    gsea_indir = os.path.abspath(f'{parent_dir}/gsea/input')
    gsea_outdir = os.path.abspath(f'{parent_dir}/gsea/output')
    
    

    def __init__(self, matrix_df, meta_df, clusters_df, error_label):
        self.matrix_df = matrix_df
        self.meta_df = meta_df
        self.clusters_df = clusters_df
        self.error_label = error_label

    def gmt_maker(self): #creates a gmt file needed for the test
        output_gmt_path = os.path.join(self.gsea_indir, 'clusters.gmt')
        with open(output_gmt_path, "w") as gmt:
            for index, row in self.clusters_df.iterrows():
                cl_name = f"cluster_{index}_{row['direction']}"
                geneline = ",".join(row["genes"].split())
                cl_desc = geneline
                clusterlist = [cl_name, cl_desc]
                samples = row["samples"].split()
                clusterlist.extend(samples)
                cluster = "\t".join(clusterlist)
                # print(repr(cluster))
                gmt.write(cluster + "\n")
        print(f"Gene sets saved to {output_gmt_path}")
        return output_gmt_path
    
    def prerank(self, output_gmt_path, metacol, samps):
        resdict_preranked = {}
        resdict_contrast = {}
        os.makedirs('gsea/input/rnks', exist_ok=True)
        rnk_outdir = os.path.abspath("gsea/input/rnks")
        rnk_df = pd.DataFrame({'Sample': self.meta_df.index, metacol : self.meta_df[metacol]})
        rnk_df = rnk_df.applymap(lambda x: None if isinstance(x, str) and self.error_label in x else x) #ERSAN Check
        rnk_cont_df = rnk_df.loc[~rnk_df.index.isin(samps), [metacol]].copy()
        rnk_path = f'{rnk_outdir}/{metacol}.rnk'
        rnk_df.to_csv(rnk_path, sep='\t', header=False, index=False)
        os.makedirs(f"GSEA/output/{metacol}__prerank", exist_ok=True)
        prerank_dirs_path = os.path.abspath(f"gsea/output/{metacol}__prerank")
        try:
            print(f'Preranking: {metacol} - Prerank')
            preranked_res = gp.prerank(rnk=rnk_path,
                threads=4, 
                gene_sets=output_gmt_path, 
                permutation_num=1000,
                min_size=1,
                outdir=None,
                no_plot=True)
            print(f'Preranked: {metacol}')
        except Exception as e:
            print(f'Error processing a column {metacol}: {e}')
        resDF_preranked = pd.DataFrame(preranked_res.res2d[['Term', 'FDR q-val']]).set_index('Term')
        try:
            print(f'Preranking: {metacol} - Contrast')
            contrast_res = gp.prerank(rnk=rnk_cont_df,
                threads=4, 
                gene_sets=output_gmt_path, 
                permutation_num=1000,
                min_size=1,
                outdir=None,
                no_plot=True)
            print(f'Preranked: {metacol}')
        except Exception as e:
            print(f'Error processing a column {metacol}: {e}')
        resDF_contrast = pd.DataFrame(contrast_res.res2d[['Term', 'FDR q-val']]).set_index('Term')
        
        return resDF_contrast, resDF_preranked

