#!/usr/bin/env python3
from sklearn.utils import resample
import pandas as pd
import numpy as np
import random
class BTSing:

    p_value_threshold = 0.05
    draws = 10


    def __init__(self, matrix_df, meta_df, cluster_df, num_metavar_list, cat_metavar_list):
        self.matrix_df = matrix_df
        self.meta_df = meta_df
        self.cluster_df = cluster_df
        self.num_metavar_list = num_metavar_list
        self.cat_metavar_list = cat_metavar_list
        # to make a dict with the cluster id as key and the sample indexes as values
        self.samples_in_cluster_dict = cluster_df['sample_indexes'].to_dict()       

        # to convert sample_indexes from string to list
        for cluster_id, sample_ids in self.samples_in_cluster_dict.items():
            self.samples_in_cluster_dict[cluster_id] = [int(sample) for sample in sample_ids.split()]




    def categorical(self):
        def metadata_per_cluster(meta_data_list, samples_in_cluster_dict) -> dict:
            metadata_per_cluster_dict = {}
            
            for cluster_id, sample_ids in samples_in_cluster_dict.items():
                values = []
                for index in sample_ids: 
                    values.append(meta_data_list[index])
                
                metadata_per_cluster_dict[cluster_id] = values
            
            return metadata_per_cluster_dict
        
        def bootstrapping(cluster_distribution_dict, draws, meta_data_list, sample_meta_data) -> dict:
            permutation_dict = {}
            length_of_cluster = len(sample_meta_data)
            
            for attribute in cluster_distribution_dict.keys():
                permutation_list = []
                
                for draw in range(draws):
                    random_draws = random.sample(meta_data_list, k=length_of_cluster)
                    attribute_counter = random_draws.count(attribute)
                    permutation_list.append(attribute_counter)
                    permutation_list.sort(reverse=True)
                    permutation_dict[attribute] = permutation_list
            
            return permutation_dict      
        # to get the distribution of sample meta data for a given cluster 
        def cluster_distribution(sample_meta_data) -> dict:
            cluster_occurrence_dict = {}
            
            for element in sample_meta_data: 
                if element in cluster_occurrence_dict:
                    cluster_occurrence_dict[element] += 1
                else:
                    cluster_occurrence_dict[element] = 1

            return cluster_occurrence_dict

        # to determine the corresponding P Values with the given cluster distribution and the bootsrapping results
        def calculate_p_value(cluster_distribution_dict, permutation_dict, draws) -> dict:
            if cluster_distribution_dict.keys() != permutation_dict.keys():
                print("keys in cluster_distribution_dict and permutation_dict don't match!")
            
            draw_freq_list = []
            p_values_dict = {}
            
            for attribute, occurrence in cluster_distribution_dict.items():
                draw_freq_list = permutation_dict[attribute]
                counter = 0
                for element in draw_freq_list:
                    if element < occurrence:
                        break
                    else:
                        counter += 1
                        
                p_values_dict[attribute] = counter/draws
            
            return p_values_dict

        # main: to harness all loops through the given categorical attributes (defined within the first cell)
        def main():
            for category in self.cat_metavar_list:
                # print(category)
                # print("------------------------------------------------------------------\n")
                meta_data_list = self.meta_df[category].tolist()
                metadata_per_cluster_dict = metadata_per_cluster(meta_data_list, self.samples_in_cluster_dict)
                
                for cluster_id, sample_meta_data in metadata_per_cluster_dict.items():
                    try:
                        cluster_distribution_dict = cluster_distribution(sample_meta_data)
                        #print(cluster_distribution_dict)
                        permutation_dict = bootstrapping(cluster_distribution_dict, self.draws, meta_data_list, sample_meta_data)
                        #print(permutation_dict)
                        p_values_dict = calculate_p_value(cluster_distribution_dict, permutation_dict, self.draws)
                        #print(p_values_dict)
                        
                        for attribute, p_value in p_values_dict.items():
                            if p_value <= self.p_value_threshold:
                               # print("Cluster ID: "+str(cluster_id)+"|", "Attribute: "+str(attribute)+"|", "P-Value: "+str(p_value)+";")
                                return p_values_dict
                    
                    except:
                        print("Exception occurred during handling of Cluster ID: "+str(cluster_id))
                        
               # print("\n") 
        p_values_dict = main()
        return p_values_dict


    def numerical(self):
        def metadata_per_cluster(meta_data_list, samples_in_cluster_dict) -> dict:
            metadata_per_cluster_dict = {}
            
            for cluster_id, sample_ids in samples_in_cluster_dict.items():
                values = []
                for index in sample_ids: 
                    values.append(meta_data_list[index])
                metadata_per_cluster_dict[cluster_id] = values
            
            return metadata_per_cluster_dict
        # to determine an aproximative distribution from the given metadata
        def bootstrapping(draws, meta_data_list, sample_meta_data):
            length_of_cluster = len(sample_meta_data)
            random_draws_mean_list = []
            
            for draw in range(draws):
                random_draws = random.sample(meta_data_list, k=length_of_cluster)
                random_draws_mean = sum(random_draws)/length_of_cluster
                random_draws_mean_list.append(random_draws_mean)
                random_draws_mean_list.sort(reverse=True)
            
            return random_draws_mean_list
        # to get the distribution of sample meta data for a given cluster 
        def cluster_distribution(sample_meta_data):
            cluster_occurrance_mean = sum(sample_meta_data)/len(sample_meta_data)
            return cluster_occurrance_mean
        # to determine the corresponding P Values with the given cluster distribution and the bootsrapping results
        def calculate_p_value(cluster_distribution_mean, permutation_list, draws) -> dict:
            counter = 0
            
            for element in permutation_list:
                if element < cluster_distribution_mean:
                    break
                else:
                    counter += 1

            p_value = counter/draws

            return p_value
        # main: to harness all loops through the given numerical attributes (defined within the first cell)
        def main():
            p_value_dict = {}
            for category in self.num_metavar_list:
                # print(category)
                # print("------------------------------------------------------------------\n")
                meta_data_list = self.meta_df[category].tolist()
                
                metadata_per_cluster_dict = metadata_per_cluster(meta_data_list, self.samples_in_cluster_dict)
                p_value_dict = {}
                for cluster_id, sample_meta_data in metadata_per_cluster_dict.items():
                    try:
                        cluster_distribution_mean = cluster_distribution(sample_meta_data)
                        #print(cluster_distribution_mean)
                        permutation_list = bootstrapping(self.draws, meta_data_list, sample_meta_data)
                        #print(permutation_list)
                        p_value = calculate_p_value(cluster_distribution_mean, permutation_list, self.draws)
                        #print(p_value)
                        
                        if p_value <= self.p_value_threshold:
                            # print("Cluster ID: "+str(cluster_id)+"|", "P-Value: "+str(p_value)+";")
                            p_value_dict[cluster_id] = p_value
                    except:
                        print("Exception occurred during handling of Cluster ID: "+str(cluster_id))
                        
                # print("\n")
            return p_value_dict
        p_values_dict = main()
        return p_values_dict





    