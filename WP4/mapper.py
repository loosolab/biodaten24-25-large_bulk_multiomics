#!/usr/bin/env python3
import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from pathlib import Path
import os

class HeatmapGen:
    parent_dir = Path(__file__).parent
    vis_outdir = os.path.abspath(f'{parent_dir}/results/heatmaps')
    def __init__(self, output_df, meta_df):
        self.output_df = output_df
        self.meta_df = meta_df 

    


    def vis_dfCreator(self, cluster: str):
        vis_df = self.output_df.drop(labels=['cluster gene ids', 'cluster genes', 'variable factors'], axis = 1)
        # print(vis_df.head())
        vis_df[['cluster', 'metacol']] = vis_df['test'].str.split(' - ', expand=True)
        vis_df.drop(columns=['test'], inplace=True)
        # vis_df_numeric = vis_df[['gsea preranked pval', 'gsea contrast pval', 'bootstrapping pval', 'ora pval']]
        vis_df['ora pval'] = pd.to_numeric(vis_df['ora pval'], errors='coerce')
        ctc = ['gsea preranked pval', 'gsea contrast pval', 'bootstrapping pval', 'ora pval']
        vis_df[ctc] = vis_df[ctc].replace('-', np.nan)
        vis_df.set_index("metacol", inplace=True)
        vis_df = vis_df[vis_df["cluster"] == cluster]

        vis_df = vis_df[ctc].apply(pd.to_numeric, errors='coerce')
        plt.figure(figsize=(12, 8))
        plt.title(cluster, fontsize=14)

        sns.heatmap(
            vis_df, 
            annot=True, 
            cmap="coolwarm", 
            fmt=".2f", 
            linewidths=0.5,
            vmax=1,
            vmin=0
            )

        plt.savefig(f"{self.vis_outdir}/{cluster}.png", dpi=300, bbox_inches="tight")
        # plt.show()
