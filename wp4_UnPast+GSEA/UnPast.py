#!/usr/bin/env python3
import os
import subprocess
from pathlib import Path

class UnPast:
    parent_dir = Path(__file__).parent
    def __init__(self, matrix_df, matrix_path, meta_df, meta_path):
        self.matrix_df = matrix_df
        self.matrix_path = matrix_path
        self.meta_df = meta_df
        self.meta_path = meta_path


    def run_unpast_docker(self, outdir, basename):
        command = [
            "docker", "run",
            "-v", f"{os.getcwd}:/user_data",
            "freddsle/unpast:latest",
            "--exprs", self.matrix_path,
            "--out_dir", outdir,
            "--basename", basename
        ]
        
        try:
            result = subprocess.run(command, check=True, text=True, capture_output=True)
            print(result.stdout)
            unpast_outfile = UnPast.unpast_outreadfiles(outdir,basename)
            return unpast_outfile
        except subprocess.CalledProcessError as e:
            print(e.stderr)

    def unpast_outreadfiles(self, directory, basename):
        listdir = os.listdir(directory)
        hitfilelist = [filename for filename in listdir if basename in filename and ".biclusters." in filename]
        if len(hitfilelist) == 0:
            print("No files were found in a given directory with a given name")
            return False
        elif len(hitfilelist) == 1:
            selectedfile = hitfilelist[0]
            print(f"Selected file: {selectedfile}")
            return selectedfile
        else:
            print("The following files fitting the given arguments were found:")
            for file in hitfilelist:
                print(file)
            while True:
                try: 
                    user_inp_selectfile = int(input("Please enter the file you want to choose: (1 for the first file, 2 for the second)  "))
                    selectedfile = hitfilelist[user_inp_selectfile-1]
                    print(f"Selected file: {selectedfile}")
                    return selectedfile
                except:
                    print("Please enter a number of the file you want to choose.")