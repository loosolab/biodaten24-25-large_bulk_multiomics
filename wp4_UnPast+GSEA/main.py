#!/usr/bin/env python3
import subprocess
import os
import pandas as pd
import argparse



def main():
    #ArgParser
    parser = argparse.ArgumentParser(description='Positional arguments: matrix path and meta path. Flags: -u for unpast, -p for preranked, -h for help.\n Usage:\b Unpast => python main.py [path/to/matrixfile.txt] path/to/metafile.txt -u [unpast/output/directory]\n Preranked => python main.py path/to/matrixfile.txt path/to/metafile.txt -p [path/to/clusterfile.txt]\n Both => python main.py [path/to/matrixfile.txt] [path/to/metafile.txt] -u [unpast/output/directory] [unpast_output_basename] -p [path/to/clusterfile.txt]')
    parser.add_argument("matrix", help='Path to the matrix file')
    parser.add_argument("meta", help='Path to the meta file Comes second after ')
    parser.add_argument("-u", nargs=2, type=str, help='Output directory and basename of the unpast output file. python main.py -u [unpast/output/directory] [unpast_output_basename] [path/to/matrixfile.txt] [path/to/metafile.txt]', default=None)
    parser.add_argument("-p", nargs=1, type=str, help='Path to the cluster file. python main.py -p [path/to/clusterfile.txt] [path/to/matrixfile.txt] [path/to/metafile.txt]', default=None,)
    # parser.add("")
    args = parser.parse_args()

    if args.u == None and args.p != None:
        prerank_run = True
        unpast_run = False
    elif args.u != None and args.p == None:
        prerank_run = False
        unpast_run = True
    elif args.u != None and args.p != None:
        prerank_run = True
        unpast_run = True
    else:
        print('Please select an analyse to run')
    

    # Initial Input
    matrix_path = os.path.abspath(args.matrix)
    meta_path = os.path.abspath(args.meta)
    print(matrix_path,meta_path)
    meta_df = pd.read_csv(meta_path, sep='\t', header=0, index_col=0, engine='python')
    matrix_df = pd.read_csv(matrix_path, sep='\t', header=0, index_col=0, engine='python')
    symbol = "|"
    matrix_df.index = matrix_df.index.str.split(symbol, n=1).str[-1]
    # print(matrix_df.head())
    # print(meta_df.head())
    
    # UnPast
    if unpast_run == True:
        try:
            unpast_outdir = os.path.abspath(args.u[0])
            unpast_basename = os.path.abspath(args.u[1])
        except Exception as e:
            print(e)
        os.makedirs(unpast_outdir, exist_ok=True)
        import UnPast as un
        unpast = un.UnPast(
            matrix_df=matrix_df,
            meta_df=meta_df,
            matrix_path=matrix_path,
            meta_path=meta_path)
        unpast_outfile = unpast.run_unpast_docker(unpast_outdir, unpast_basename)
        unpast_outpath="/".join([unpast_outdir, unpast_outfile])

    #Prerank
    if prerank_run == True:
        if unpast_run == True:
            clusters_df = pd.read_csv(unpast_outpath, sep="\t", header=0, index_col=0)
        else:
            unpast_outpath = os.path.abspath(args.p[0])
            clusters_df = pd.read_csv(unpast_outpath, sep="\t", header=0, index_col=0)
        import GSEA as gs
        gsea = gs.GSEA(
            matrix_df=matrix_df,
            meta_df=meta_df,
            clusters_df=clusters_df)
        gmt_file_path = gsea.gmt_maker()
        gsea.prerank(gmt_file_path)



# print(clusters_df.head())





if __name__ == "__main__":
    main()