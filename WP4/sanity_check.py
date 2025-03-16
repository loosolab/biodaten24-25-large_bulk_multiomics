#!/usr/bin/env python3

# """
# This script performs a sanity check on a metadata file by detecting and marking outliers and invalid values.
# The script uses the following steps:
# 1. Load the metadata file from *meta_data_source*.
# 2. Iterate over all attributes and check for outliers or invalid values.
# 3. For categorical attributes, replace invalid values with "ERRSAN" (categorical attributes must be defined under *categorical_attributes* before execution).
# 4. For binary attributes, replace invalid values with "ERRSAN" (binary attributes must be defined under *binary_attributes* before execution).
# 5. For numerical attributes, detect outliers using Isolation Forest and replace them with "ERRSAN" (all values that are not categorical or binary are considered numerical!). 
# 6. Save the updated metadata file under *output_path*.

# Usage:
# Run the script with the appropriate paths to the metadata file and output file. The script will process the metadata and save the updated file with outliers and invalid values marked as "ERRSAN".
# """

import pandas as pd
import numpy as np
from sklearn.ensemble import IsolationForest
from pathlib import Path

# Path to the meta data as well as the output file and definition of the categorical and binary attributes


# Definition of the outlier detection parameters
class SanCheck:
    paparent_dir = Path(__file__).parent

    categorical_attributes = {"gender":["F", "M"], "condition": ["healthy","pah","ph-lung"]}
    binary_attributes = ["Paradoxe_Septumbewegung", "Perikarderguss"]

    def __init__(self, determination_threshold, number_of_estimators, meta_data_df):
        self.determination_threshold = determination_threshold   # 4.4 is based on various testing and should be adjusted according to the data
        self.number_of_estimators = number_of_estimators     # Number of estimators for the Isolation Forest
        self.meta_data_df = meta_data_df
        



    def outlier_detection(self, column):
        # Convert data to a NumPy array
        column_array = np.array(column).reshape(-1, 1)

        # Isolation Forest with automatic adjustment of anomaly detection
        iforest = IsolationForest(n_estimators=self.number_of_estimators, contamination='auto', random_state=42)
        iforest.fit(column_array)

        # Calculate score (the lower, the more likely an outlier)
        scores = iforest.decision_function(column_array)

        # Adaptive threshold for outlier determination
        threshold = np.mean(scores) - self.determination_threshold * np.std(scores)

        # Mark as outlier
        outliers = (scores < threshold).astype(int)
        outlier_indices = np.where(outliers == 1)[0]
        
        # Mark the outliers in the column
        for idx in outlier_indices:
            column[idx] = f"ERRSAN[{column[idx]}]"
        
        return column

    def SanityCheck(self):
        # Load the meta data file and assign the column names
        meta_data_df = self.meta_data_df
        column_names = meta_data_df.columns.tolist()[1:]

        # Iterate over all attributes and check for outliers or invalid values
        for n, attribute in enumerate(column_names):
            # Get the list of values for the current attribute and print the progress
            meta_data_list = meta_data_df[attribute].tolist()
            print(f"Processing attribute: {n+1}/{len(column_names)} ({attribute})")
            
            # Check if the attribute is categorical and replace invalid values with ERRSAN
            if attribute in self.categorical_attributes.keys():
                for i, value in enumerate(meta_data_list):
                    if value not in self.categorical_attributes[attribute] or value == "NaN":
                        meta_data_list[i] = f"ERRSAN[{value}]"
                meta_data_df[attribute] = meta_data_list
            
            # Check if the attribute is binary and replace invalid values with ERRSAN
            elif attribute in self.binary_attributes:
                for i, value in enumerate(meta_data_list):
                    if pd.isna(value):
                        meta_data_list[i] = f"ERRSAN[{value}]"
                    elif value not in [0.0, 1.0]:
                        meta_data_list[i] = f"ERRSAN[{value}]"
                    else:
                        meta_data_list[i] = float(value)
                meta_data_df[attribute] = meta_data_list
            
            # When attribute is numerical, check for outliers and replace them with ERRSAN
            else:
                errsan_dict = {}
                non_errsan_list = []
                
                # Check for NaN and invalid values and replace them with ERRSAN
                for i, value in enumerate(meta_data_list):
                    if pd.isna(value) or (isinstance(value, str) and not value.startswith("ERRSAN")):
                        meta_data_list[i] = f"ERRSAN[{value}]"
                    
                # Sort the values into previous assigned ERRSAN values and numerical values
                for i, value in enumerate(meta_data_list):
                    if isinstance(value, str) and value.startswith("ERRSAN"):
                        errsan_dict[i] = value
                    else:
                        non_errsan_list.append(float(value))
                
                # Check for outliers in the non_errsan_list and replace them with ERRSAN, merge with errsan_dict
                outlier_list = self.outlier_detection(column=non_errsan_list)
                for idx, value in errsan_dict.items():
                    outlier_list.insert(idx, value)
            
                # Update the meta data file with the new values
                meta_data_df[attribute] = outlier_list

        # Save the updated meta data file
        #print(meta_data_df)
        return meta_data_df
