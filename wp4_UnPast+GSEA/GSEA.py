#!/usr/bin/env python3
import os
import pandas as pd
import gseapy as gp
from pathlib import Path
class GSEA:
    parent_dir = Path(__file__).parent
    def __init__(self, matrix_df, meta_df, clusters_df):
        self.matrix_df = matrix_df
        self.meta_df = meta_df
        self.clusters_df = clusters_df
    def gmt_maker(self):
        output_gmt_path = os.path.join(self.parent_dir, "GSEA/clusters.gmt")
        with open(output_gmt_path, "w") as gmt:
            # symbol = "|"
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
    def prerank(self, output_gmt_path):
        meta_int_columns = []
        for column in self.meta_df.columns:
            if pd.to_numeric(self.meta_df[column].dropna(), errors='coerce').notna().all(): # weird way to check for numeric columns
                meta_int_columns.append(column)
        # print(meta_int_columns)
        os.makedirs('GSEA/output/rnks', exist_ok=True)
        rnk_outdir = os.path.abspath("GSEA/output/rnks")
        for col in meta_int_columns:
            rnk_df = pd.DataFrame({'Sample': self.meta_df.index, col : self.meta_df[col]})
            # print(rnk_df.head())
            rnk_path = f'{rnk_outdir}/{col.replace('/',' ')}.rnk'
            rnk_df.to_csv(rnk_path, sep='\t', header=False, index=False)
            os.makedirs(f"GSEA/output/{col.replace('/',' ')}_prerank", exist_ok=True)
            prerank_dirs_path = os.path.abspath(f"GSEA/output/{col.replace('/',' ')}_prerank")
            try:
                print(f'Preranking: {col}')
                gp.prerank(rnk=rnk_path,
                    threads=4, 
                    gene_sets=output_gmt_path, 
                    permutation_num=1000, 
                    outdir=prerank_dirs_path)
                print(f'Preranked: {col}')
            except Exception as e:
                print(f'Error processing a column {col}: {e}')





