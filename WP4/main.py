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
    numMetaVar_configList = config['numerical_columns'].split(',')
    numMetaVar_configList = [e.strip().lower() for e in numMetaVar_configList]
    catMetaVar_configList = config['categorical_columns'].split(',')
    catMetaVar_configList = [e.strip().lower() for e in catMetaVar_configList]
    error_label = config['error_label']
    results_path = f'{parent_dir}/results/{config['results_name']}'

    # Initial Input, matrix and meta dataframes
    matrix_path = os.path.abspath(f'{parent_dir}/init_input/{matrixfile_name}')
    meta_path = os.path.abspath(f'{parent_dir}/init_input/{metafile_name}')
    meta_df = pd.read_csv(meta_path, sep='\t', header=0, index_col=0, engine='python')
    matrix_df = pd.read_csv(matrix_path, sep='\t', header=0, index_col=0, engine='python')

    # Sanity Check -> should produce the same meta_df but edited with ERSANs
    import sanity_check as snch
    sanCheck = snch.SanCheck(
        determination_threshold = 4.4,   # 4.4 is based on various testing and should be adjusted according to the data
        number_of_estimators = 1000,
        meta_data_df = meta_df
    )
    meta_df = sanCheck.SanityCheck()
    symbol = "|"
    matrix_df.index = matrix_df.index.str.split(symbol, n=1).str[-1]          # separating and removing ENS.... gene identifiers
    meta_df.columns = meta_df.columns.str.replace('/', '&')
    print(testsdict)

    num_metavar_list = []
    cat_metavar_list = []

    for column in meta_df.columns:
        if column.lower().strip() in numMetaVar_configList:
            num_metavar_list.append(column)
        elif column.lower().strip() in catMetaVar_configList:
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
    # ANALYSIS of PVALs
    # GSEA GMT Maker
    if testsdict['gsea'] == True:
        import GSEA as gs
        gsea = gs.GSEA(
            matrix_df=matrix_df,
            meta_df=meta_df,
            clusters_df=clusters_df,
            error_label = error_label)
        gmt_file_path = gsea.gmt_maker()
        # gsea.prerank(gmt_file_path)
    # gsea_outdir = f'{parent_dir}/gsea/output'
    # contrast_df = gsea.contrast(gmt_file_path)

    # Bootstrapping
    if testsdict['bootstrapping'] == True:
        import Bootstrapping as bts
        bts = bts.BTSing()
    
    #ORA
    if testsdict['ora'] == True:
        import ORA as ora
        ora = ora.ORA()

    # Final results (creation of the .csv table)
    # iterationg through: 1. metadata table, 2. cluster table from unpast
    columns_list = ['test', 'cluster genes', 'cluster gene ids', 'variable factors', 'gsea preranked pval', 'gsea contrast pval' , 'bootstrapping pval', 'ora pval']
    # columns_list = [ele for ele in columns_list if testsdict.get(ele.split(' ')[0], None) != False]
    output_df = pd.DataFrame(columns=columns_list)
    result_list = []
    # Column Iteration
    for column in meta_df.columns:
        print(f'Iteration is now on column: {column}')
        if column in num_metavar_list:
            num_meta_var = column             #num_meta_var
            column_iscat = False                              #if the column is categorical
        elif column in cat_metavar_list:
            cat_meta_var = column              #cat_meta_var
            column_iscat = True
        else:
            column_iscat = None
            continue
        # GSEA
        if testsdict['gsea'] == True and column_iscat == False:
            # prerank
            gsea_pval_df_list = gsea.prerank(output_gmt_path=gmt_file_path, 
                                            metacol=column, 
                                            samps = cl_samples)
            gseaPvalPreranked_df = gsea_pval_df_list[0]
            gseaPvalContrast_df = gsea_pval_df_list[1]
        # Cluster Iteration
        for cl_index, cl_row in clusters_df.iterrows():
            cl_name = f"cluster_{cl_index}_{cl_row['direction']}"    #cl name
            print(f'Iteration is now on cluster: {cl_name}')
            cl_samples = cl_row['samples'].split()
            cl_gene_num = cl_row['n_genes']                  #n genes
            cl_gene_ids = cl_row['genes']                    #gene names
            # I. categorical
            # 1. BTS
            if column_iscat == True:
                unique_vars_list = meta_df[column].unique().tolist()
                unique_vars_count = len(unique_vars_list)   
                for unvar in unique_vars_list:
                    # 1. BTS 
                    if testsdict['bootstrapping'] == True:
                        cat_bts_pval_dict = bts.categorical(meta_column=meta_df[column], cl_samples=cl_samples)
                        bts_pval = cat_bts_pval_dict[unvar]
                    # 2. ORA
                    if testsdict['ora'] == True:
                        cat_ora_pval_dict = ora.categorical(meta_column=meta_df[column], cl_samples=cl_samples)
                        ora_pval = cat_ora_pval_dict[unvar]
                    
                    # creating and appending the resultlist as a line in final output table
                    # values of the 'test' column differ depending on whether the metavar is categotical or numerical
                    # for nums the result line will be appended once for each metadata of each cluster
                    # for cats the result line will be appended for each instance of metadata column (male, female etc.)

                    # Resultlist cat
                    result_list.append({
                        'test': f"{cl_name} - {cat_meta_var}/{unvar}".replace('&','/'),
                        'cluster genes': cl_gene_num,
                        'cluster gene ids': cl_gene_ids,
                        'variable factors': unique_vars_count,
                        'gsea preranked pval': '-',
                        'gsea contrast pval': '-',
                        'bootstrapping pval': bts_pval if 'bts_pval' in locals() else None,
                        'ora pval': ora_pval if 'ora_pval' in locals() else None
                    })
                    bts_pval = None
                    ora_pval = None
                    print('Cat added to the list')
            # II. numerical
            elif column_iscat == False:
                # 1. GSEA
                if testsdict['gsea'] == True:
                    gseaPreranked_pval = gseaPvalPreranked_df.loc[cl_name, 'FDR q-val'] if cl_name in gseaPvalPreranked_df.index else None
                    gseaContrast_pval = gseaPvalContrast_df.loc[cl_name, 'FDR q-val'] if cl_name in gseaPvalContrast_df.index else None
                # 2. BTS
                if testsdict['bootstrapping'] == True:
                    bts_pval = bts.numerical(meta_column=meta_df[column], cl_samples=cl_samples)
            # Resultlist num
            if column in num_metavar_list:
                result_list.append({
                    'test': f"{cl_name} - {num_meta_var}".replace('&','/'),
                    'cluster genes': cl_gene_num,
                    'cluster gene ids': cl_gene_ids,
                    'variable factors': 'NA',
                    'gsea preranked pval': gseaPreranked_pval if 'gseaPreranked_pval' in locals() else 'Err172',
                    'gsea contrast pval': gseaContrast_pval if 'gseaContrast_pval' in locals() else 'Err173',
                    'bootstrapping pval': bts_pval if 'bts_pval' in locals() else '?Err174',
                    'ora pval': '-'
                })
                gseaPreranked_pval = None
                gseaContrast_pval = None
                bts_pval = None
                ora_pval = None
                print('Num added to the list')
    # print(f'Final list: {result_list}')
    output_df = pd.concat([output_df, pd.DataFrame(result_list)], ignore_index=True)
    output_df.to_csv(results_path, index=False)
    print(f"Final results saved as {results_path}")

            


                        


if __name__ == "__main__":
    main()