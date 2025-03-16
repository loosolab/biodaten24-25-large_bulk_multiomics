#!/usr/bin/env python3
import os
import pandas as pd
import gseapy as gp
from pathlib import Path
import numpy as np

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
        self.output_gmt_path = self.gmt_maker()
    
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

    def prerank(self, metacol):
        output_gmt_path = self.output_gmt_path
        rnk_df = pd.DataFrame({'Sample': self.meta_df.index, metacol : self.meta_df[metacol]})
        rnk_df = rnk_df.applymap(lambda x: None if isinstance(x, str) and self.error_label in x else x) #ERSAN Check
        rnk_df.set_index('Sample', inplace=True)
        try:
            print(f'Preranking: {metacol} - Prerank')
            print(f'Prerank RNK: ')
            print(rnk_df)
            print('Gene sets: ')
            print(output_gmt_path)
            preranked_res = gp.prerank(rnk=rnk_df,
                threads=4, 
                gene_sets=output_gmt_path, 
                permutation_num=1000,
                min_size=1,
                outdir=None,
                no_plot=True,
                verbose=True)
            print(f'Preranked: {metacol}')
        except Exception as e:
            print(f'Error processing a column {metacol}: {e}')
        resDF_preranked = pd.DataFrame(preranked_res.res2d[['Term', 'FDR q-val']]).set_index('Term')
        print('Prerank Result:')
        print(resDF_preranked)
        return resDF_preranked

    def contrast(self, metacol, samps, cluster):
        output_gmt_path = self.output_gmt_path     
        rnk_df = pd.DataFrame({'Sample': self.meta_df.index, metacol : self.meta_df[metacol]})
        rnk_df = rnk_df.applymap(lambda x: None if isinstance(x, str) and self.error_label in x else x) #ERSAN Check
        rnk_df.set_index('Sample', inplace=True)
        rnk_cont_df = rnk_df.loc[~rnk_df.index.isin(samps), [metacol]].copy()
        dg = "DUMMY_GENE"     # To avoid the cluster being filtered out from geneset
        rnk_cont_df.loc[dg] = {'Sample' : dg, metacol: round(rnk_df[metacol].mean(), 1)}
        with open(output_gmt_path, "r") as infile:
            for line in infile:
                if cluster in line:
                    parts = line.strip().split("\t")
                    parts.append(dg)
                    genline = "\t".join(parts[0:1])
                    gendict = {genline: parts[2:]}
                    print(f'Gendict: {gendict}')
                    try:
                        print(f'Preranking: {metacol} - Contrast, cluster: {cluster}')
                        print(f'Samples {samps}')
                        print(f'Contrast RNK: ')
                        print(rnk_cont_df)
                        contrast_res = gp.prerank(rnk=rnk_cont_df,
                            threads=4,
                            gene_sets=gendict,
                            permutation_num=1000,
                            min_size=1,
                            outdir=None,
                            no_plot=True,
                            verbose=True)
                        print(f'Contrasted: {metacol}')
                    except Exception as e:
                        print(f'Error processing a column {metacol}: {e}')
        resDF_contrast = pd.DataFrame(contrast_res.res2d[['Term', 'FDR q-val']]).set_index('Term')
        print('Contrast:')
        print(resDF_contrast)
        return resDF_contrast

