#!/usr/bin/env python3
from pathlib import Path
import os
import pandas as pd
import argparse
import json



def main():
    parent_dir = Path(__file__).parent
    #ArgParser
    parser = argparse.ArgumentParser(description='Matrix and meta file should be in "init_input" folder. Copy the name of your meta and matrix file in the corresponding fields in config.json')
    parser.add_argument("-n", nargs=1, type=str, help='Type the names of analyses to skip them. Seprate with "_". For example: unpast_gsea_bootstrapping_ora will skip all possible analyses', default=None)    
    parser.add_argument("matrixname", type=str, help='The name of the matrix file. The file should be in "init_input" folder', default=None)
    parser.add_argument("metaname", type=str, help='The name of the meta file. The file should be in "init_input" folder', default=None)
    

    args = parser.parse_args()
    testsdict = {"unpast":True,"gsea":True,'bootstrapping':True,'ora':True}
    blocklist = str(args.n[0]).lower().strip().split("_")
    for test in blocklist:
        if testsdict[test] == True:       # check the undesired tests via -n. 
            testsdict[test] = False    

    # config json
    config_path = f'{parent_dir}/config.json'
    with open(config_path, "r") as file:
        config = json.load(file)

    matrixfile_name = args.matrixname
    metafile_name = args.metaname
    unpast_basename = config['unpast_basename']
    cluserfile_name = config['clusterfile_name']
    results_path = f'{parent_dir}/results/{config['results_name']}'

    # Initial Input, matrix and meta dataframes
    matrix_path = os.path.abspath(f'{parent_dir}/init_input/{matrixfile_name}')
    meta_path = os.path.abspath(f'{parent_dir}/init_input/{metafile_name}')
    meta_df = pd.read_csv(meta_path, sep='\t', header=0, index_col=0, engine='python')
    matrix_df = pd.read_csv(matrix_path, sep='\t', header=0, index_col=0, engine='python')
    symbol = "|"
    matrix_df.index = matrix_df.index.str.split(symbol, n=1).str[-1]          # separating and removing ENS.... gene identifiers
    meta_df.columns = meta_df.columns.str.replace('/', '&')
    print(testsdict)

    num_metavar_list = []
    cat_metavar_list = []
    for column in meta_df.columns:
        if pd.to_numeric(meta_df[column].dropna(), errors='coerce').notna().all():
            num_metavar_list.append(column)
        elif pd.to_numeric(meta_df[column].dropna(), errors='coerce').notna().all() == False and not 'id' in column:
            cat_metavar_list.append(column)
        else:
            continue

    # UnPast
    if testsdict['unpast'] == True:
        import UnPast as un
        unpast = un.UnPast(
            matrix_df=matrix_df,
            meta_df=meta_df,
            matrix_path=matrix_path,
            meta_path=meta_path)
        unpast_outfile = unpast.run_unpast_docker(unpast_basename)
        unpast_outpath="/".join([unpast.unpast_outdir, unpast_outfile])

    # Clusterfile DataFrame creation 
    if testsdict['unpast'] == True:
        clusterpath = unpast_outpath
    else:
        clusterpath = f'{parent_dir}/clusters/{cluserfile_name}'
    clusters_df = pd.read_csv(clusterpath, sep="\t", header=0, index_col=0)

    # GSEA Prerank
    if testsdict['gsea'] == True:
        import GSEA as gs
        gsea = gs.GSEA(
            matrix_df=matrix_df,
            meta_df=meta_df,
            clusters_df=clusters_df)
        gmt_file_path = gsea.gmt_maker()
        gsea.prerank(gmt_file_path)
    gsea_outdir = f'{parent_dir}/gsea/output'

    # Bootstrapping
    if testsdict['bootstrapping'] == True:
        import Bootstrapping as bts
        bts = bts.BTSing()
    
    #ORA
    if testsdict['ora'] == True:
        import ORA as ora
        ora = ora.ORA()

    # Final results (creation of the .csv table)
    # iterationg through: 1. cluster table from unpast, 2. metadata table
    columns_list = ['test', 'cluster genes', 'cluster gene ids', 'variable factors', 'gsea pval', 'bootstrapping pval', 'ora pval']
    columns_list = [ele for ele in columns_list if testsdict.get(ele.split(' ')[0], None) != False]
    output_df = pd.DataFrame(columns=columns_list)
    result_list = []
    for cl_index, cl_row in clusters_df.iterrows():
        cl_name = f"cluster_{cl_index}_{cl_row['direction']}"    #cl name
        print(f'Cl name: {cl_name}')
        cl_samples = cl_row['samples'].split()
        cl_gene_num = cl_row['n_genes']                  #n genes
        cl_gene_ids = cl_row['genes']                    #gene names
        for column in meta_df.columns:
            print(f'column: {column}')
            if column in num_metavar_list:
                num_meta_var = column             #num_meta_var
                column_iscat = False                              #if the column is categorical
            elif column in cat_metavar_list:
                cat_meta_var = column              #cat_meta_var
                column_iscat = True
            else:
                column_iscat = None
                continue
            # GSEA pval
            if testsdict['gsea'] == True and column_iscat == False:
                for dirnames in os.listdir(gsea_outdir):  
                    if num_meta_var.lower() == dirnames.split('__')[0].lower() and Path(os.path.join(gsea_outdir, dirnames)).is_dir() and '__prerank' in dirnames:
                        gsea_outfile = os.path.join(gsea_outdir, dirnames, "gseapy.gene_set.prerank.report.csv")
                    else: 
                        continue
                    if os.path.exists(gsea_outfile):
                        gsea_outdf = pd.read_csv(gsea_outfile)
                        for index,row in gsea_outdf.iterrows():
                            if cl_name == row['Term'] and not row.empty:
                                gsea_pval = row["FDR q-val"]  # FDR q-val
            else:
                gsea_pval = None
            # categorical
            # bootstrapping pval
            if column_iscat == True:
                unique_vars_list = meta_df[column].unique().tolist()
                unique_vars_count = len(unique_vars_list)   
                for unvar in unique_vars_list:
                    if testsdict['bootstrapping'] == True:
                        cat_bts_pval_dict = bts.categorical(meta_column=meta_df[column], cl_samples=cl_samples)
                        bts_pval = cat_bts_pval_dict[unvar]
                    if testsdict['ora'] == True:
                        cat_ora_pval_dict = ora.categorical(meta_column=meta_df[column], cl_samples=cl_samples)
                        ora_pval = cat_ora_pval_dict[unvar]
                    result_list.append({
                        'test': f"{cl_name} - {cat_meta_var}/{unvar}".replace('&','/'),
                        'cluster genes': cl_gene_num,
                        'cluster gene ids': cl_gene_ids,
                        'variable factors': unique_vars_count,
                        'gsea pval': '-',
                        'bootstrapping pval': bts_pval if 'bts_pval' in locals() else None,
                        'ora pval': ora_pval if 'ora_pval' in locals() else None
                    })
                    bts_pval = None
                    ora_pval = None
                    print('Cat added to the list')
            # numerical
            elif testsdict['bootstrapping'] == True and column_iscat == False:
                bts_pval = bts.numerical(meta_column=meta_df[column], cl_samples=cl_samples)

            # creating and appending the resultlist as a line in final output table
            # values of the 'test' column differ depending on whether the metavar is categotical or numerical
            # for nums the result line will be appended once for each metadata of each cluster
            # for cats the result line will be appended for each instance of metadata column (male, female etc.)
            if column in num_metavar_list:
                result_list.append({
                    'test': f"{cl_name} - {num_meta_var}".replace('&','/'),
                    'cluster genes': cl_gene_num,
                    'cluster gene ids': cl_gene_ids,
                    'variable factors': 'NA',
                    'gsea pval': gsea_pval if 'gsea_pval' in locals() else None,
                    'bootstrapping pval': bts_pval if 'bts_pval' in locals() else None,
                    'ora pval': ora_pval if 'ora_pval' in locals() else None
                })
                gsea_pval = None
                bts_pval = None
                ora_pval = None
                print('Num added to the list')
    # print(f'Final list: {result_list}')
    output_df = pd.concat([output_df, pd.DataFrame(result_list)], ignore_index=True)
    output_df.to_csv(results_path, index=False)
    print(f"Final results saved as {results_path}")

            


                        


if __name__ == "__main__":
    main()