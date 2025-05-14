import re
import pandas as pd
import numpy as np

from quatradis.gene.gene import Gene
from quatradis.embl.reader import EMBLReader
from quatradis.tisp.parser import PlotParser
import time
import math


class GeneAnnotator:
    def __init__(self,old_algorithm, annotation_file,blocks,plotfiles_all=None,forward_count_condition=None,reverse_count_condition=None,combined_count_condition=None,forward_count_control=None,reverse_count_control=None,combined_compare_csv=None,forward_compare_csv=None,reverse_compare_csv=None,**kwargs):
        
        self.annotation_file = annotation_file
        self.embl_reader = EMBLReader(self.annotation_file)
        if old_algorithm:
            self.blocks = self.sort_blocks_by_start_coord(blocks)
            self.knockout_proportion_start = 0.5
            self.increased_expression_proportion_end = 0.3
            self.features = self.embl_reader.features
        else:
            self.condition_files, self.control_files = self.get_condition_control_separately(plotfiles_all)
            self.forward_count_df= self.read_combine_count_files(forward_count_condition)
            self.reverse_count_df=self.read_combine_count_files(reverse_count_condition)
            self.forward_count_control_df= self.read_combine_count_files(forward_count_control)
            self.reverse_count_control_df=self.read_combine_count_files(reverse_count_control)
            self.combined_count_condition_df_rep1= self.read_tsv(combined_count_condition[0])
            self.combined_compare_csv= combined_compare_csv
            self.forward_compare_csv= forward_compare_csv
            self.reverse_compare_csv= reverse_compare_csv
            self.embl_reader = EMBLReader(self.annotation_file)
            self.features = self.embl_reader.genes_to_features
            self.merged_forward_reverse_compare_csv=None
            self.result_df=None
            self.blue_red_logfc_diff_threshold = kwargs.get("blue_red_logfc_diff_threshold", None)
            self.distance_threshold = kwargs.get("distance_threshold", None)
            self.insertion_count_max_threshold = kwargs.get("insertion_count_max_threshold", None)
            self.insertion_signal_similarity_avg_threshold = kwargs.get("insertion_signal_similarity_avg_threshold", None)
            self.log_fc_threshold = kwargs.get("log_fc_threshold", None)
            self.overlap_threshold = kwargs.get("overlap_threshold", None)
            self.q_val_threshold = kwargs.get("q_val",None)
            self.inactivation_fraction_threshold = kwargs.get("inactivation_fraction_threshold",None)
            self.condition_plots= self.read_plot_files(self.condition_files)
            self.control_plots=self.read_plot_files(self.control_files)
            total_read_condition_plots= [plot.total_reads for plot in self.condition_plots]
            print("total_read_condition_plots",total_read_condition_plots)
            average_total_condition_reads= sum(total_read_condition_plots)/len(total_read_condition_plots)
            print("average_total_condition_reads",average_total_condition_reads)
            reference_total_reads= 33346524.0
            self.insertion_count_sum_threshold = math.ceil(kwargs.get("insertion_count_sum_threshold", None)*len(self.condition_files) * (average_total_condition_reads/reference_total_reads))
            print("self.insertion_count_sum_threshold",self.insertion_count_sum_threshold)
       
    def read_plot_files(self, files):
        parsed_plots_list=[]
        for file in files:
            parsed_plots_list.append(PlotParser(file))
        return parsed_plots_list
    
    def read_tsv(self,file):
        return pd.read_csv(file,sep="\t")

    def get_condition_control_separately(self,plotfiles_all):
        condition_plot_files=[]
        control_plot_files=[]
        num_files = len(plotfiles_all)

        if num_files % 2 != 0:
            raise ValueError("The total number of plotfiles (control + conditions) must be equal.")

        half_index = num_files // 2

        for i in range(half_index):
            condition_file = plotfiles_all[i]
            control_file = plotfiles_all[i + half_index]

            condition_plot_files.append(condition_file)
            control_plot_files.append(control_file)

        return condition_plot_files, control_plot_files

    def sort_blocks_by_start_coord(self, blocks):
        sorted_blocks = sorted((b for b in blocks), key=lambda x: x.start)
        return sorted_blocks
        
    def get_inactivation_fraction(self, start,end,strand, condition_plots):
        """
        Computes the inactivation fraction for a gene based on the insertions in both the forward and reverse 
        strands across different conditions, normalized by the gene length.

        The inactivation fraction is calculated by determining the first index in the gene's region 
        where all values across conditions are positive (for both strands). If no such index exists, 
        it returns the first positive value across all conditions. The result is then normalized by the gene's length.

        Args:
        - f (object): A feature object for the gene (not directly used in the current method but could 
                    be useful for future extension).
        - gene (Gene): A gene object that contains `start` and `end` attributes representing the 
                    genomic coordinates of the gene.
        - condition_plots (list): A list of objects (e.g., PlotParser) representing the data for different 
                                conditions, each with `forward` and `reverse` attributes, which store 
                                insertion data for the forward and reverse strands.

        Returns:
        - float: The inactivation fraction for the gene, normalized by its length. The result is a value 
                between 0 and 1, where a result of -1 indicates that no valid index was found.

        Example:
        - gene = Gene(start=100, end=200)
        condition_plots = [PlotParser(file1), PlotParser(file2), ...]
        inactivation_fraction = get_inactivation_fraction(f, gene, condition_plots)
        This will compute and return the inactivation fraction based on the gene's insertions across conditions.
        """
        gene_length = end - start

        # Slice the arrays for both forward and reverse strands
        forward_gene_inserts_array = [arr.forward[start:end] for arr in condition_plots]
        reverse_gene_inserts_array = [arr.reverse[start:end] for arr in condition_plots]

        if strand == -1:
            forward_gene_inserts_array = [arr[::-1] for arr in forward_gene_inserts_array]
            reverse_gene_inserts_array = [arr[::-1] for arr in reverse_gene_inserts_array]
  
        # Combine into matrices
        forward_sliced_matrix = np.vstack(forward_gene_inserts_array)
        reverse_sliced_matrix = np.vstack(reverse_gene_inserts_array)

        # Helper function to calculate the result index for a given matrix
        def find_result_index(sliced_matrix):
            # Find the first positive value index across all arrays
            first_positive_indices = [
                np.where(arr > 5)[0][0] if np.any(arr > 5) else float("inf") for arr in sliced_matrix
            ]
            return min(first_positive_indices)

        # Find result indices for both forward and reverse strands
        forward_result_index = find_result_index(forward_sliced_matrix)
        reverse_result_index = find_result_index(reverse_sliced_matrix)

        # Get the minimum of the two result indices
        result_index = min(forward_result_index, reverse_result_index)
            # If the result is infinity, return -1
        # if result_index == float("inf"):
        #     result_index = -1

        # Normalize the result by gene length
        return 1 - (result_index) / gene_length

      
    
    def get_genes_with_up_or_down_regulation(self):
        """
        This function retrieves the gene names that have either upregulation or downregulation based on 
        the 'Category2' and 'Category3' columns of the result dataframe. 

        The function performs the following steps:
        - Filters the result dataframe to include only the rows where 'Category2' or 'Category3' is not null.
        - Extracts the unique gene names from the filtered dataframe.

        Returns:
        - unique_gene_names: A list of unique gene names that are upregulated or downregulated based on 
        the conditions set in the 'Category2' and 'Category3' columns.
        """
        filtered_df = self.result_df[self.result_df['Category2'].notnull()| self.result_df['Category3'].notnull() ]
        unique_gene_names = filtered_df['Gene'].unique()
        return unique_gene_names
    
    def check_insertion_signal_pattern_shape(self, signals, distance_threshold=0.6):
        """
        Check if multiple signals share a similar pattern based on Shape-Based Distance (SBD),
        handling zero-variance signals.

        Parameters:
            signals (list of numpy.ndarray): List of 1D signals to compare.
            distance_threshold (float): Threshold for similarity (between 0 and 1).

        Returns:
            bool: True if all valid signals have pairwise distances below the threshold, False otherwise.
            dict: Pairwise shape-based distance values for analysis.
        """
        num_signals = len(signals)
        distances = {}  # To store pairwise shape-based distances

        # Compute first-order differences and normalize each signal
        processed_signals = []
        valid_indices = []
        for i, signal in enumerate(signals):
            signal = np.array(signal)
            if len(signal) < 2:
                print(f"Signal {i+1} is invalid (too short). Skipping.")
                processed_signals.append(None)
                continue

            diff_signal = np.diff(signal)  # First-order difference
            norm = np.linalg.norm(diff_signal)

            if norm == 0:  # Handle zero-variance signal
                print(f"Signal {i+1} is invalid (zero variance). Skipping.")
                processed_signals.append(None)
            else:
                processed_signals.append(diff_signal / norm)
                valid_indices.append(i)
        
        # Compute pairwise shape-based distances
        for i in range(num_signals):
            for j in range(i + 1, num_signals):
                key = f"Signal_{i+1}_vs_Signal_{j+1}"
                
                # Skip if either signal is invalid
                if processed_signals[i] is None or processed_signals[j] is None:
                    distances[key] = None  # Mark as invalid comparison
                    continue
                
                # Ensure both signals have the same length (trim longer one)
                min_len = min(len(processed_signals[i]), len(processed_signals[j]))
                signal1, signal2 = processed_signals[i][:min_len], processed_signals[j][:min_len]

                # Compute Shape-Based Distance (SBD)
                dist = np.linalg.norm(signal1 - signal2)
                distances[key] = dist
                print(f"SBD between Signal {i+1} and Signal {j+1}: {dist:.4f}")

        # Check if all valid pairwise distances are below the threshold
        all_similar = all(dist is not None and dist <= distance_threshold for dist in distances.values())

        return all_similar, distances


    def check_average_threshold(self,lists):
        averages = []
        max_avg = 0

        # Compute non-zero averages and track max in a single pass
        for lst in lists:
            non_zero = [x for x in lst if x != 0]
            if non_zero:
                avg = sum(non_zero) / len(non_zero)
                averages.append(avg)
                max_avg = max(max_avg, avg)


        # Check if all averages are within X% of max_avg
        threshold = self.insertion_signal_similarity_avg_threshold* max_avg
        # print("averages:",averages)
        return all(avg >= threshold for avg in averages)


    def insertion_signal_similarity_check(self, check,geneName=None):
        """
        This function checks the similarity of insertion signals for a given gene and its corresponding 
        3' or 5' prime regions in condition files. It compares the insertion signals based on the gene's strand and location 
        across multiple condition files.

        The function performs the following steps:
        - Initializes the condition plots for each condition file.
        - Determines the appropriate insertion signal array type (either forward or reverse) based on the 
        gene's strand and the provided check (either "5prime" or "3prime").
        - Retrieves the insertion signals for the gene from the specified region.
        - Calls the check_insertion_signal_pattern method to compare the insertion signals and calculate 
        pairwise correlations.

        Args:
        - f: Gene object containing gene details, including location and strand.
        - check: A string indicating whether to check the "5prime" or "3prime" region.

        Returns:
        - all_similar: A boolean indicating if the insertion signal pattern is consistent across conditions.
        
        Raises:
        - ValueError: If an invalid combination of strand and check is provided.
        
        Note:
        - This function is used to assess insertion signal consistency before further categorization of gene activity.
        """
        ## Takes time - recode it.
        # condition_plots=[]
        # # if not hasattr(self, "condition_plots"):
        # for file in self.condition_files:
        #     condition_plots.append(PlotParser(file))
        condition_plots= self.condition_plots
        start = int(self.combined_count_condition_df_rep1.loc[self.combined_count_condition_df_rep1["gene_name"] == geneName, "start"].iloc[0])
        end = int(self.combined_count_condition_df_rep1.loc[self.combined_count_condition_df_rep1["gene_name"] == geneName, "end"].iloc[0])
        strand = int(self.combined_count_condition_df_rep1.loc[self.combined_count_condition_df_rep1["gene_name"] == geneName, "strand"].iloc[0])
        mapping = {
        (1, "5prime"): "forward",
        (-1, "5prime"): "reverse",
        (1, "3prime"): "reverse",
        (-1, "3prime"): "forward",
        }
    
        # Get the appropriate attribute based on the mapping
        array_type = mapping.get((strand, check))
        if array_type:
            gene_inserts_array = [
                getattr(arr, array_type)[start:end] for arr in condition_plots
            ]
        else:
            raise ValueError("Invalid combination of strand and check")

        # all_similar_correlation, pairwise_correlations= self.check_insertion_signal_pattern(gene_inserts_array,-self.insertion_lag,self.insertion_lag)
        all_similar_distances, pairwise_distances= self.check_insertion_signal_pattern_shape(gene_inserts_array,self.distance_threshold)
        all_similar_averages= self.check_average_threshold(gene_inserts_array)
        all_similar= all_similar_distances and all_similar_averages
        # if geneName=="marA__5prime":
        #     print(geneName)
        #     print("pairwise_correlations",pairwise_correlations)
        #     print("all_similar_averages",all_similar_averages)
        #     print("all_similar",all_similar)
        #     print("****************")
        print(geneName)
        print("pairwise_distances",pairwise_distances)
        print("all_similar_distances",all_similar_distances)
        print("all_similar_averages",all_similar_averages)
        print("all_similar",all_similar)
        print("****************")
        return all_similar
    
    
    def categorize_activitation_explanation(self,geneName):
        # f = self.features[geneName]
        strand = int(self.combined_count_condition_df_rep1.loc[self.combined_count_condition_df_rep1["gene_name"] == geneName, "strand"].iloc[0])
        five_prime_name = f"{geneName}__5prime"
        three_prime_name = f"{geneName}__3prime"
        if strand == 1:
            logfc_key_3_same, logfc_key_3_opp = "logFC_reverse", "logFC_forward"
            logfc_key_5_same, logfc_key_5_opp = "logFC_forward", "logFC_reverse"
        else:
            logfc_key_3_same, logfc_key_3_opp = "logFC_forward", "logFC_reverse"
            logfc_key_5_same, logfc_key_5_opp = "logFC_reverse", "logFC_forward"

        #Check for 5prime/UpRegulation
        
        prime_3 = self.merged_forward_reverse_compare_csv[self.merged_forward_reverse_compare_csv["gene_name"] == three_prime_name]
        prime_5 = self.merged_forward_reverse_compare_csv[self.merged_forward_reverse_compare_csv["gene_name"] == five_prime_name]
        if not prime_5.empty:
            all_similar= self.insertion_signal_similarity_check("5prime",five_prime_name)

            if all_similar:
                logfc_same = prime_5[logfc_key_5_same].values
                logfc_opp= prime_5[logfc_key_5_opp].values
                if not pd.isna(logfc_same) and not pd.isna(logfc_opp):
                    prime5_overlaps= self.embl_reader.find_overlaps(five_prime_name)
                    if len(prime5_overlaps)>0:
                        for gene , overlap_pct in prime5_overlaps:
                            
                            # if overlap_pct<self.overlap_threshold:
                            #     if abs(logfc_same.item()-logfc_opp.item())<self.blue_red_logfc_diff_threshold:
                            #         if five_prime_name=="marA__5prime":
                            #             print("overlap_pct",overlap_pct)
                            #             print("logfc_same.item()-logfc_opp.item():",logfc_same.item()-logfc_opp.item())
                            #         self.result_df.loc[self.result_df['Gene'] == geneName, 'Category2'] = None

                            # else:
                            if overlap_pct>self.overlap_threshold:
                                gene_knockout= self.result_df[self.result_df["Gene"]==gene]["Category1"].values
                                if abs(logfc_same.item()-logfc_opp.item()<self.blue_red_logfc_diff_threshold):
                                    new_value=f"upregulated due to knockout of gene {gene}"
                                    if gene_knockout.size > 0 and gene_knockout[0] is not None:
                                        if gene_knockout[0]=="knockout":
                                            self.result_df.loc[self.result_df['Gene'] == geneName, 'Category4'] = (
                                                self.result_df.loc[self.result_df['Gene'] == geneName, 'Category4']
                                                .apply(lambda x: new_value if pd.isna(x) else f"{x} | {new_value}")
                                            )
                                elif logfc_same.item()> logfc_opp.item()+ self.blue_red_logfc_diff_threshold:
                                    new_value=f"Knockout due to upregulation of gene {geneName}"
                                    # if not pd.isna(gene_knockout):
                                    if gene_knockout.size > 0 and gene_knockout[0] is not None:
                                        if gene_knockout[0]=="knockout":
                                        # if gene_knockout.item()=="knockout":
                                            self.result_df.loc[self.result_df['Gene'] == gene, 'Category4'] = (
                                                self.result_df.loc[self.result_df['Gene'] == gene, 'Category4']
                                                .apply(lambda x: new_value if pd.isna(x) else f"{x} | {new_value}")
                                            )
            else:
                # if five_prime_name=="marA__5prime":
                #     print(five_prime_name,"-",all_similar)

                self.result_df.loc[self.result_df['Gene'] == geneName, 'Category2'] = None
        
        if not prime_3.empty:
            all_similar= self.insertion_signal_similarity_check("3prime",three_prime_name)
            if all_similar:
        
                logfc_same = prime_3[logfc_key_3_same].values
                logfc_opp= prime_3[logfc_key_3_opp].values
                if not pd.isna(logfc_same) and not pd.isna(logfc_opp):
                    prime3_overlaps= self.embl_reader.find_overlaps(three_prime_name)
                    if len(prime3_overlaps)>0:
                        for gene , overlap_pct in prime3_overlaps:
                            # if overlap_pct<self.overlap_threshold:
                            #     if abs(logfc_same.item()-logfc_opp.item()<self.blue_red_logfc_diff_threshold):
                            #         self.result_df.loc[self.result_df['Gene'] == geneName, 'Category3'] = None
                            # else:
                            if overlap_pct>self.overlap_threshold:
                                gene_knockout= self.result_df[self.result_df["Gene"]==gene]["Category1"].values
                                if abs(logfc_same.item()-logfc_opp.item())<self.blue_red_logfc_diff_threshold:
                                    new_value=f"downregulated due to knockout of gene {gene}"
                                    # if not pd.isna(gene_knockout):
                                    if gene_knockout.size > 0 and gene_knockout[0] is not None:
                                        # if gene_knockout.item()=="knockout":
                                        if gene_knockout[0] == "knockout":
                                            self.result_df.loc[self.result_df['Gene'] == geneName, 'Category4'] = (
                                                self.result_df.loc[self.result_df['Gene'] == geneName, 'Category4']
                                                .apply(lambda x: new_value if pd.isna(x) else f"{x} | {new_value}")
                                        )
                                elif logfc_same.item()> logfc_opp.item() + self.blue_red_logfc_diff_threshold:
                                    new_value=f"Knockout due to downregulation of gene {geneName}"
                                    if gene_knockout.size > 0 and gene_knockout[0] is not None:
                                        if gene_knockout[0] == "knockout":
                                            self.result_df.loc[self.result_df['Gene'] == gene, 'Category4'] = (
                                                self.result_df.loc[self.result_df['Gene'] == gene, 'Category4']
                                                .apply(lambda x: new_value if pd.isna(x) else f"{x} | {new_value}")
                                            ) 
            else:
                
                self.result_df.loc[self.result_df['Gene'] == geneName, 'Category3'] = None
                    
    
    def get_unique_genes_that_have_significant_3_5_prime_logfc(self):
        """
        Cleans the gene names by removing '3prime', '5prime', and any optional '__', then filters 
        for genes with significant logFC values (qval<5%). Returns a list of unique 
        gene names that meet these criteria.

        Returns:
        --------
        list
            A list of unique gene names with significant logFC values (qval<5%) after 
            cleaning '3prime' and '5prime' markers.
        """
        df = pd.DataFrame()
        # Cleaning the gene names by removing '3prime' and '5prime' along with optional '__'
        df['gene_name_cleaned'] = self.merged_forward_reverse_compare_csv['gene_name'].str.replace(r'(__)?(3prime|5prime)', '', regex=True)
        unique_gene_names = df['gene_name_cleaned'].dropna().unique()
        unique_gene_names = list(unique_gene_names)
        return unique_gene_names
    
    def read_combine_count_files(self,file_paths):
        count_dfs = [pd.read_csv(file_path, sep='\t') for file_path in file_paths]
        final_df = (
                    pd.concat(count_dfs)            # Combine all DataFrames
                    .groupby('gene_name', as_index=False)  # Group by 'gene_name'
                    .agg({'read_count': 'sum', 'ins_index': 'sum'})     # Sum 'read_count' for each 'gene_name'
                )
        return final_df
    
    def insertion_threshold_check(self,result_df):
        """
        Filters and updates rows in the given DataFrame based on insertion count thresholds.

        This function applies filtering conditions based on gene strand direction and 
        regulation categories. It performs the following operations:

        1. Merges forward and reverse read count data to obtain insertion counts.
        2. Determines whether genes meet the insertion count threshold based on their strand 
        and regulatory category (upregulated/downregulated).
        3. Flags rows where insertion counts fall below the threshold.
        4. Updates or removes flagged rows:
        - If a row is flagged for both upregulated (Category2) and downregulated (Category3), it is removed.
        - If only Category2 is flagged, its value along with 'LogFC(5_Prime)' is set to NaN.
        - If only Category3 is flagged, its value along with 'LogFC(3_Prime)' is set to NaN.

        Parameters:
        - result_df (pd.DataFrame): The input DataFrame containing gene, strand, and 
        regulatory category information.

        Returns:
        - pd.DataFrame: A cleaned and updated DataFrame with filtered entries based on 
        insertion count thresholds.
        """
        #Condition analysis count files
        forward_count_df = self.forward_count_df.rename(columns={'read_count': 'forward_count', 'ins_index': 'forward_ins_index'})
        reverse_count_df = self.reverse_count_df.rename(columns={'read_count': 'reverse_count', 'ins_index': 'reverse_ins_index'})
        merged_count_df= pd.merge(forward_count_df, reverse_count_df, on='gene_name', how='inner')

        #Control analysis count files
        forward_count_control_df = self.forward_count_control_df.rename(columns={'read_count': 'forward_count', 'ins_index': 'forward_ins_index'})
        reverse_count_control_df = self.reverse_count_control_df.rename(columns={'read_count': 'reverse_count', 'ins_index': 'reverse_ins_index'})
        merged_count_control_df= pd.merge(forward_count_control_df, reverse_count_control_df, on='gene_name', how='inner')

        # dictionary for fast lookups - condition
        forward_count_dict = merged_count_df.set_index('gene_name')['forward_count'].to_dict()
        reverse_count_dict = merged_count_df.set_index('gene_name')['reverse_count'].to_dict()
        forward_ins_index_dict = merged_count_df.set_index('gene_name')['forward_ins_index'].to_dict()
        reverse_ins_index_dict = merged_count_df.set_index('gene_name')['reverse_ins_index'].to_dict()

        # dictionary for fast lookups - control
        forward_count_control_dict = merged_count_control_df.set_index('gene_name')['forward_count'].to_dict()
        reverse_count_control_dict = merged_count_control_df.set_index('gene_name')['reverse_count'].to_dict()
        forward_ins_index_control_dict = merged_count_control_df.set_index('gene_name')['forward_ins_index'].to_dict()
        reverse_ins_index_control_dict = merged_count_control_df.set_index('gene_name')['reverse_ins_index'].to_dict()
        
        # Ensure Category2 and Category3 are treated as strings and handle NaN values
        result_df['Category1'] = result_df['Category1'].astype(str).fillna('')
        result_df['Category2'] = result_df['Category2'].astype(str).fillna('')
        result_df['Category3'] = result_df['Category3'].astype(str).fillna('')

        # Initialize flags
        flag_cat2 = []
        flag_cat3 = []
        flag_remove_knockout = []
        flag_remove_protection=[]
        insertion_count_upregulated = []
        insertion_count_downregulated = []
        insertion_index_upregulated = []
        insertion_index_downregulated = []
        insertion_count_max_upregulated=[]
        insertion_count_max_downregulated=[]
        condition_plots=[]
        control_plots=[]
        

        # for file in self.condition_files:
        #     condition_plots.append(PlotParser(file))
        condition_plots= self.condition_plots

    
        # for file in self.control_files:
        #     control_plots.append(PlotParser(file))
        control_plots= self.control_plots

        test_results={}
        for _, row in result_df.iterrows():
            gene = row['Gene']
            strand = row['Strand']
            category1_knockout = row['Category1'] == 'knockout'
            category1_protection = row['Category1'] == 'protection'
            category2_up = row['Category2'] == 'upregulated'
            category3_down = row['Category3'] == 'downregulated'
            
            gene_5prime = f"{gene}__5prime"
            gene_3prime = f"{gene}__3prime"
            
            cat2_flag = False
            cat3_flag = False
            remove_knockout_flag = False
            remove_protection_flag= False
            
            if category2_up:
                if strand == 1:
                    cat2_count_flag = forward_count_dict.get(gene_5prime, 0) < self.insertion_count_sum_threshold
                    # f = self.features[gene_5prime]
                    start = int(self.combined_count_condition_df_rep1.loc[self.combined_count_condition_df_rep1["gene_name"] == gene_5prime, "start"].iloc[0])
                    end = int(self.combined_count_condition_df_rep1.loc[self.combined_count_condition_df_rep1["gene_name"] == gene_5prime, "end"].iloc[0])
                    # strand = int(self.combined_count_condition_df_rep1.loc[self.combined_count_condition_df_rep1["gene_name"] == gene_5prime, "strand"].iloc[0])
                    gene_inserts_array = [getattr(arr, "forward")[start:end] for arr in condition_plots]
                    cat2_count_mean= np.mean([np.max(np.array(arr))  for arr in gene_inserts_array]) 
                    cat2_max_flag= cat2_count_mean < self.insertion_count_max_threshold
                    cat2_flag= cat2_count_flag or cat2_max_flag
                    insertion_count_upregulated.append(forward_count_dict.get(gene_5prime, 0))
                    insertion_index_upregulated.append(forward_ins_index_dict.get(gene_5prime, 0))
                    insertion_count_max_upregulated.append(cat2_count_mean)
                    test_results[gene_5prime]=cat2_flag

                elif strand == -1:
                    cat2_count_flag = reverse_count_dict.get(gene_5prime, 0) < self.insertion_count_sum_threshold
                    # f = self.features[gene_5prime]
                    start = int(self.combined_count_condition_df_rep1.loc[self.combined_count_condition_df_rep1["gene_name"] == gene_5prime, "start"].iloc[0])
                    end = int(self.combined_count_condition_df_rep1.loc[self.combined_count_condition_df_rep1["gene_name"] == gene_5prime, "end"].iloc[0])

                    gene_inserts_array = [getattr(arr, "reverse")[start:end] for arr in condition_plots]
                    cat2_count_mean= np.mean([np.max(np.array(arr))  for arr in gene_inserts_array]) 
                    cat2_max_flag= cat2_count_mean < self.insertion_count_max_threshold
                    cat2_flag= cat2_count_flag or cat2_max_flag
                    insertion_count_upregulated.append(reverse_count_dict.get(gene_5prime, 0))
                    insertion_index_upregulated.append(reverse_ins_index_dict.get(gene_5prime, 0))
                    insertion_count_max_upregulated.append(cat2_count_mean)
                    test_results[gene_5prime]=cat2_flag
        
            else:
                insertion_count_upregulated.append(np.nan)
                insertion_index_upregulated.append(np.nan)
                insertion_count_max_upregulated.append(np.nan)
            if category3_down:
                if strand == 1:
                    cat3_count_flag = reverse_count_dict.get(gene_3prime, 0) < self.insertion_count_sum_threshold
                    # f = self.features[gene_3prime]
                    start = int(self.combined_count_condition_df_rep1.loc[self.combined_count_condition_df_rep1["gene_name"] == gene_3prime, "start"].iloc[0])
                    end = int(self.combined_count_condition_df_rep1.loc[self.combined_count_condition_df_rep1["gene_name"] == gene_3prime, "end"].iloc[0])
                    gene_inserts_array = [getattr(arr, "reverse")[start:end] for arr in condition_plots]
                    cat3_count_mean= np.mean([np.max(np.array(arr)) for arr in gene_inserts_array]) 
                    cat3_max_flag= cat3_count_mean < self.insertion_count_max_threshold
                    cat3_flag= cat3_count_flag or cat3_max_flag
                    insertion_count_max_downregulated.append(cat3_count_mean)
                    insertion_count_downregulated.append(reverse_count_dict.get(gene_3prime, 0))
                    insertion_index_downregulated.append(reverse_ins_index_dict.get(gene_3prime, 0))
                    test_results[gene_3prime]=cat3_flag
                elif strand == -1:
                    cat3_count_flag = forward_count_dict.get(gene_3prime, 0) < self.insertion_count_sum_threshold
                    # f = self.features[gene_3prime]
                    start = int(self.combined_count_condition_df_rep1.loc[self.combined_count_condition_df_rep1["gene_name"] == gene_3prime, "start"].iloc[0])
                    end = int(self.combined_count_condition_df_rep1.loc[self.combined_count_condition_df_rep1["gene_name"] == gene_3prime, "end"].iloc[0])
                    gene_inserts_array = [getattr(arr, "forward")[start:end] for arr in condition_plots]
                    cat3_count_mean= np.mean([np.max(np.array(arr)) for arr in gene_inserts_array]) 
                    cat3_max_flag= cat3_count_mean < self.insertion_count_max_threshold
                    cat3_flag= cat3_count_flag or cat3_max_flag
                    insertion_count_max_downregulated.append(cat3_count_mean)
                    insertion_count_downregulated.append(forward_count_dict.get(gene_3prime, 0))
                    insertion_index_downregulated.append(forward_ins_index_dict.get(gene_3prime, 0))
                    test_results[gene_3prime]=cat3_flag
            else:
                insertion_count_downregulated.append(np.nan)
                insertion_index_downregulated.append(np.nan)
                insertion_count_max_downregulated.append(np.nan)

            if category1_knockout:
                total_count = forward_count_dict.get(gene, 0) + reverse_count_dict.get(gene, 0)
                remove_knockout_sum_flag = total_count < 1.25*self.insertion_count_sum_threshold
                # f = self.features[gene]
                start = int(self.combined_count_condition_df_rep1.loc[self.combined_count_condition_df_rep1["gene_name"] == gene, "start"].iloc[0])
                end = int(self.combined_count_condition_df_rep1.loc[self.combined_count_condition_df_rep1["gene_name"] == gene, "end"].iloc[0])
                gene_inserts_array = [getattr(arr, "combined")[start:end] for arr in condition_plots]
                cat1_count_mean= np.mean([np.max(np.array(arr)) for arr in gene_inserts_array]) 
                remove_knockout_max_flag= cat1_count_mean < self.insertion_count_max_threshold
                remove_knockout_flag= remove_knockout_sum_flag or remove_knockout_max_flag
                test_results[gene]=remove_knockout_flag
            if category1_protection:
                total_count = forward_count_control_dict.get(gene, 0) + reverse_count_control_dict.get(gene, 0)
                remove_protection_sum_flag = total_count < 1.25*self.insertion_count_sum_threshold
                # f= self.features[gene]
                start = int(self.combined_count_condition_df_rep1.loc[self.combined_count_condition_df_rep1["gene_name"] == gene, "start"].iloc[0])
                end = int(self.combined_count_condition_df_rep1.loc[self.combined_count_condition_df_rep1["gene_name"] == gene, "end"].iloc[0])
                gene_inserts_array = [getattr(arr, "combined")[start:end] for arr in control_plots]
                cat1_count_mean= np.mean([np.max(np.array(arr)) for arr in gene_inserts_array])
                remove_protection_max_flag= cat1_count_mean < self.insertion_count_max_threshold
                remove_protection_flag= remove_protection_sum_flag or remove_protection_max_flag
                test_results[gene]=remove_protection_flag

            
            flag_cat2.append(cat2_flag)
            flag_cat3.append(cat3_flag)
            flag_remove_knockout.append(remove_knockout_flag)
            flag_remove_protection.append(remove_protection_flag)

        
        result_df['flag_cat2'] = flag_cat2
        result_df['flag_cat3'] = flag_cat3
        result_df['flag_remove_knockout'] = flag_remove_knockout
        result_df['flag_remove_protection']= flag_remove_protection
        result_df['insertion_count_upregulated'] = insertion_count_upregulated
        result_df['insertion_count_downregulated'] = insertion_count_downregulated
        result_df['insertion_index_upregulated'] = insertion_index_upregulated
        result_df['insertion_index_downregulated'] = insertion_index_downregulated
        result_df['insertion_count_max_upregulated'] = insertion_count_max_upregulated
        result_df['insertion_count_max_downregulated'] = insertion_count_max_downregulated
        # print("test_results")
        # print(test_results)

        # with open("insert_count_check_data.json", "w") as f:
        #     json.dump(test_results, f, indent=4)
        
        # Drop rows where all flags are True
        result_df = result_df[~(result_df['flag_cat2'] & result_df['flag_cat3'] & result_df['flag_remove_knockout'] & result_df['flag_remove_protection'])]

        
        # Update Category2 and LogFC(5_Prime) if flag_cat2 is True and flag_cat3 is False
        mask_cat2 = (result_df['flag_cat2']) & (~result_df['flag_cat3'])
        result_df.loc[mask_cat2, ['Category2', 'LogFC(5_Prime)']] = np.nan
        
        # Update Category3 and LogFC(3_Prime) if flag_cat3 is True and flag_cat2 is False
        mask_cat3 = (result_df['flag_cat3']) & (~result_df['flag_cat2'])
        result_df.loc[mask_cat3, ['Category3', 'LogFC(3_Prime)']] = np.nan

        mask_knockout = result_df['flag_remove_knockout']
        result_df.loc[mask_knockout, ['Category1','LogFC(Gene)']] = np.nan

        mask_protection = result_df['flag_remove_protection']
        result_df.loc[mask_protection, ['Category1','LogFC(Gene)']] = np.nan
        
        # Drop flag columns before returning
        result_df.drop(columns=['flag_cat2', 'flag_cat3', 'flag_remove_knockout','flag_remove_protection'], inplace=True)
        
        return result_df.reset_index(drop=True)
    
    def signal_strength_score(self, df):
        """
        Calculate confidence scores for genes based on their regulation status.

        This function computes two confidence scores:
        1. `confidence_score_upregulated` for rows where 'Category2' is 'upregulated'.
        2. `confidence_score_downregulated` for rows where 'Category3' is 'downregulated'.

        The confidence score is calculated using the following formulas:

        - If 'Category2' == 'upregulated':
        confidence_score_upregulated = (LogFC(5_Prime) ** 0.5) * 
                                        ((insertion_index_upregulated / number_replicates) ** 0.1) * 
                                        ((insertion_count_upregulated / number_replicates) ** 0.5)

        - If 'Category3' == 'downregulated':
        confidence_score_downregulated = (LogFC(3_Prime) ** 0.5) * 
                                        ((insertion_index_downregulated / number_replicates) ** 0.1) * 
                                        ((insertion_count_downregulated / number_replicates) ** 0.5)

        Parameters:
        -----------
        df : pd.DataFrame
            Input DataFrame containing the required columns:
            - 'Category2' (for upregulation classification)
            - 'Category3' (for downregulation classification)
            - 'LogFC(5_Prime)', 'LogFC(3_Prime)'
            - 'insertion_index_upregulated', 'insertion_count_upregulated'
            - 'insertion_index_downregulated', 'insertion_count_downregulated'

        Returns:
        --------
        pd.DataFrame
            The DataFrame with two new columns:
            - 'confidence_score_upregulated': Computed confidence score for upregulated genes.
            - 'confidence_score_downregulated': Computed confidence score for downregulated genes.

        Notes:
        ------
        - The variable `number_replicates` is determined by the length of `self.condition_files`.
        - Rows that do not satisfy the conditions ('upregulated' for Category2 or 'downregulated' for Category3) 
        will have NaN values in the corresponding confidence score column.
        """

        number_replicates= len(self.condition_files)

        # Apply the formula only for 'upregulated' rows, otherwise set NaN
        df["confidence_score_upregulated"] = np.where(
            df["Category2"] == "upregulated",
            (df["LogFC(5_Prime)"] ** 0.1) *
            ((df["insertion_index_upregulated"] / number_replicates) ** 0.1) *
            ((df["insertion_count_upregulated"] / number_replicates) ** 0.5),
            np.nan  # Assign NaN for non-'upregulated' rows
        )

        df["confidence_score_downregulated"] = np.where(
            df["Category3"] == "downregulated",
            (df["LogFC(3_Prime)"] ** 0.5) *
            ((df["insertion_index_downregulated"] / number_replicates) ** 0.1) *
            ((df["insertion_count_downregulated"] / number_replicates) ** 0.5),
            np.nan  # Assign NaN for non-'upregulated' rows
        )

        return df        
    
    def categorize_knockout_protection(self, gene_combined_csv):
        """
        Categorizes genes based on their log fold change (logFC) and other criteria.

        This function processes gene-level data from the combined CSV and determines 
        a gene's category (e.g., "knockout," "protection," or "unclassified") based 
        on logFC thresholds and transcription fraction criteria. It also extracts 
        relevant gene features and metrics.

        Args:
            gene_combined_csv (pd.Series): A row of the combined CSV containing gene-specific 
                                        data such as 'gene_name', 'logFC', and 'q.value'.

        Returns:
            list: A list containing:
                - Gene name
                - Primary category
                - Secondary and tertiary categories (None for now)
                - Gene start and end positions
                - Gene strand
                - Log fold change (logFC)
                - LogFC for 3' prime and 5' prime (None for now)
        
        Notes:
            - This function is under development and aims to support confidence score generation 
            for gene categorizations.
            - Categories are determined based on logFC thresholds:
                - Positive logFC above a thershold and with inactivation percentage below the threshold: "knockout"
                - Negative logFC: "protection"
                - Otherwise: "unclassified"
            - Placeholder values are returned if the gene name matches a specific format.
        """
        gene_name = gene_combined_csv["gene_name"]

        # Skip genes with "__3prime" or "__5prime"
        if not re.search(r"^(.+)__([35])prime$", gene_name):
            # Precompute plots only once if they don't change per row
            # if not hasattr(self, "condition_plots"):
            #     self.condition_plots = [PlotParser(file) for file in self.condition_files]

            start = int(self.combined_count_condition_df_rep1.loc[self.combined_count_condition_df_rep1["gene_name"] == gene_name, "start"].iloc[0])
            end = int(self.combined_count_condition_df_rep1.loc[self.combined_count_condition_df_rep1["gene_name"] == gene_name, "end"].iloc[0])
            strand = int(self.combined_count_condition_df_rep1.loc[self.combined_count_condition_df_rep1["gene_name"] == gene_name, "strand"].iloc[0])
            logfc_gene = gene_combined_csv["logFC"]
            qval_gene = gene_combined_csv["q.value"]
            logcpm = gene_combined_csv["logCPM"]

            # Determine category
            if logfc_gene > self.log_fc_threshold:
                inactivation_fraction = self.get_inactivation_fraction(start,end,strand, self.condition_plots)
                gene_category = "fractional knockout" if inactivation_fraction <= self.inactivation_fraction_threshold else "knockout"
            elif logfc_gene < -self.log_fc_threshold:
                gene_category = "protection"
            else:
                gene_category = "unclassified"

            return [
                gene_name, gene_category, None, None, start, end, strand,
                logfc_gene, None, None, qval_gene, None, None, logcpm, None, None, None, None
            ]
        else:
            return [None] * 18

    def read_and_filter_logfc_qval(self,csv_file_path,categorization_type=None):
        """
        Reads a CSV file containing gene expression data, filters the data to retain only the rows 
        where the 'q.value' column is less than 0.05, and returns the filtered dataframe.

        The function performs the following steps:
        - Loads the CSV file into a pandas dataframe.
        - Filters the dataframe to include only rows where the 'q.value' column has a value less than 0.05 and 'LogFC' >threshold or <-threshold.

        Args:
        - csv_file_path (str): The path to the CSV file containing the gene expression data.

        Returns:
        - df (pandas.DataFrame): A filtered dataframe containing only the rows where the 'q.value' 
        column is less than 0.05.

        Note:
        - The filtering is performed to retain only the data with statistically significant results (based on a 
        q-value threshold of default value(0.05) and LogFC > threshold or less than -LogFC).
        """
        
        df= pd.read_csv(csv_file_path)

        if categorization_type== "prime_ends":
            df_filtered = df[(df["q.value"] < self.q_val_threshold)]
        else:
            df_filtered = df[(df["q.value"] < self.q_val_threshold) & ((df["logFC"] > self.log_fc_threshold) | (df["logFC"] < -self.log_fc_threshold))]

        return df_filtered

    def prime_end_expression_categorization(self, geneName):
        """
        Modifies the regulation categorization of a given gene by analyzing its 3' and 5' prime regions 
        and calculating log fold changes (logFC) to categorize the gene as upregulated or downregulated 
        based on specific thresholds.

        The function performs the following steps:
        - Determines the strand direction (forward or reverse) for the gene.
        - Fetches the log fold change values for the 3' and 5' prime regions of the gene.
        - Categorizes the gene as upregulated or downregulated based on the log fold change values and a predefined threshold.
        - Updates the `result_df` DataFrame with the regulation categories and log fold changes for the gene.

        Args:
        - geneName (str): The name of the gene for which the regulation categorization will be modified.

        Returns:
        - None: The function modifies the `result_df` DataFrame directly.

        Note:
        - The function assumes that the `result_df` DataFrame already contains gene-related data and 
        will update it with the new regulation categorization information.
        - The function checks if the gene is present in the `result_df` and either updates the existing row 
        or adds a new row if the gene is not found.

        Example:
        - geneName = "GeneA"
        After calling `prime_end_expression_categorization("GeneA")`, the function will update the 
        regulation categories (up/downregulated) and log fold changes in `result_df` for "GeneA".
        """
        start = int(self.combined_count_condition_df_rep1.loc[self.combined_count_condition_df_rep1["gene_name"] == geneName, "start"].iloc[0])
        end = int(self.combined_count_condition_df_rep1.loc[self.combined_count_condition_df_rep1["gene_name"] == geneName, "end"].iloc[0])
        strand = int(self.combined_count_condition_df_rep1.loc[self.combined_count_condition_df_rep1["gene_name"] == geneName, "strand"].iloc[0])
        # f = self.features[geneName]
        # strand = f.strand

        def calculate_logfc_and_category(prime_3, prime_5, is_forward):
            nonlocal up_regulation_category, down_regulation_category,logfc_3_prime,logfc_5_prime, Mlogfc,qval_3_prime, qval_5_prime,logCPM_3_prime, logCPM_5_prime

            if not prime_3.empty:
                prime3_logfc = prime_3[logfc_key_3].values
                if not pd.isna(prime3_logfc) and prime3_logfc > self.log_fc_threshold:
                    down_regulation_category = "downregulated"
                    logfc_3_prime = prime3_logfc.item()
                    qval_3_prime= prime_3[qval_key_3].values.item()
                    logCPM_3_prime= prime_3[logCPM_key_3].values.item()
                    # Mlogfc = (prime3_logfc * prime3_ins_idx[prime3_ins_idx['gene_name']==three_prime_name]['min_ins_index'].values[0]) / prime_3[qval_key_3].values
            if not prime_5.empty:
                prime5_logfc = prime_5[logfc_key_5].values
                if not pd.isna(prime5_logfc) and prime5_logfc > self.log_fc_threshold:
                    up_regulation_category = "upregulated"
                    logfc_5_prime = prime5_logfc.item()
                    qval_5_prime= prime_5[qval_key_5].values.item()
                    logCPM_5_prime= prime_5[logCPM_key_5].values.item()
            
        # Determine keys based on strand
        logfc_key_3 = "logFC_reverse" if strand == 1 else "logFC_forward"
        logfc_key_5 = "logFC_forward" if strand == 1 else "logFC_reverse"
        qval_key_3 = "q.value_reverse" if strand == 1 else "q.value_forward"
        qval_key_5 = "q.value_forward" if strand == 1 else "q.value_reverse"
        logCPM_key_3 = "logCPM_reverse" if strand == 1 else "logCPM_forward"
        logCPM_key_5 = "logCPM_forward" if strand == 1 else "logCPM_reverse"
        # prime3_ins_idx = self.reverse_ins_idx if strand == 1 else self.forward_ins_idx
        # prime5_ins_idx = self.forward_ins_idx if strand == 1 else self.reverse_ins_idx

        # Fetch prime data
        five_prime_name = f"{geneName}__5prime"
        three_prime_name = f"{geneName}__3prime"
        prime_3 = self.merged_forward_reverse_compare_csv[self.merged_forward_reverse_compare_csv["gene_name"] == three_prime_name]
        prime_5 = self.merged_forward_reverse_compare_csv[self.merged_forward_reverse_compare_csv["gene_name"] == five_prime_name]

        up_regulation_category, down_regulation_category,logfc_3_prime,logfc_5_prime, Mlogfc,qval_3_prime, qval_5_prime, logCPM_3_prime, logCPM_5_prime = None, None, None,None,None,None,None,None,None
        calculate_logfc_and_category(prime_3, prime_5, strand == 1)
        if not isinstance(self.result_df, pd.DataFrame):
            self.result_df = pd.DataFrame(columns=['Gene', 'Category1', 'Category2', 'Category3', 'Start', 'End', 
                                                'Strand', 'LogFC(Gene)', 'LogFC(3_Prime)', 'LogFC(5_Prime)',
                                                'Qval(Gene)', 'Qval(5_Prime)', 'Qval(3_Prime)', 
                                                'Log_CPM(Gene)', 'Log_CPM(5_Prime)', 'Log_CPM(3_Prime)'])

        if geneName in self.result_df["Gene"].values:
            self.result_df.loc[self.result_df["Gene"] == geneName, ['Category2','Category3', 'LogFC(3_Prime)', 'LogFC(5_Prime)','Qval(3_Prime)','Qval(5_Prime)','Log_CPM(3_Prime)','Log_CPM(5_Prime)']] = [
                up_regulation_category, down_regulation_category,logfc_3_prime, logfc_5_prime, qval_3_prime,qval_5_prime,logCPM_3_prime,logCPM_5_prime]
        else:
            new_row = {
                'Gene': geneName,
                'Category1':None,
                'Category2': up_regulation_category,
                'Category3':down_regulation_category,
                'Start': start,
                'End': end,
                'Strand': strand,
                'LogFC(Gene)': None,
                'LogFC(3_Prime)':logfc_3_prime,
                 'LogFC(5_Prime)':logfc_5_prime,
                 'Qval(Gene)':None,
                 'Qval(5_Prime)':qval_5_prime,
                 'Qval(3_Prime)':qval_3_prime,
                 'Log_CPM(Gene)':None,
                 'Log_CPM(5_Prime)':logCPM_5_prime,
                 'Log_CPM(3_Prime)':logCPM_3_prime
            }
            self.result_df = pd.concat([self.result_df, pd.DataFrame([new_row])], ignore_index=True)
    
    def annotate_genes_new(self):
        """
            Annotates genes based on provided data and generates a DataFrame with categorized results.

            This method performs the following steps:
            1. Reads and filters logFC (log Fold Change) values from combined, forward, and reverse compare csv files.
            2. Applies activation categorization to the combined CSV data.
            3. Merges forward and reverse compare csv files to identify unique gene names.
            4. Categorizes genes based on modified regulation and activation explanations.
            5. Cleans the resulting DataFrame, replaces missing values, and returns the annotated results.
            6. Inspect Inflation of Logfc for low insertions (Inprogress). 

            Returns:
                pandas.DataFrame: A DataFrame containing annotated genes with the following columns:
                    - 'Gene'
                    - 'Category1'       (Knockout | Protected | Unclassified categorization of gene)
                    - 'Category2'       (Upregulated | None categorization of gene)
                    - 'Category3'       (Downregulated | None categorization of gene)
                    - 'Start'           (Starting index of gene)
                    - 'End'             (Ending index of gene)
                    - 'Strand'          (Location/Position of gene)
                    - 'LogFC(Gene)'     (LogFC of gene with qval<5% if gene is categorized as knockout or protected else None)
                    - 'LogFC(3_Prime)'  (LogFC of 3 prime with qval<5% if gene is categorized as downregualted else None)
                    - 'LogFC(5_Prime)'  (LogFC of 5 prime with qval<5% if gene is categorized as downregualted else None)
                    - 'Category4' (added during processing explaining if the upregualtion is due to knockout of adjacent gene or vice versa.)

            Example:
                result_df = gene_annotator.annotate_genes()
        """
        time_qval_filter= time.time()
        # filter genes based on logfc significance (qval<5%) for combined compare insertions
        self.combined_compare_csv= self.read_and_filter_logfc_qval(self.combined_compare_csv)
        # print("Time to filter Qvals",time.time()-time_qval_filter)
        # perform gene knockout, protection catgeorization
        time_category1_categorization= time.time()
        print("Result_DF before knockout-protection",self.result_df)
        print("combined_compare_csv", self.combined_compare_csv)
        knockout_protection_result_list=[]
        for idx, row in self.combined_compare_csv.iterrows():
            result = self.categorize_knockout_protection(row)
            knockout_protection_result_list.append(result)
        # self.result_df= self.combined_compare_csv.apply(self.categorize_knockout_protection, axis=1, result_type='expand')
      
        print("Error Here")
        self.result_df = pd.DataFrame(knockout_protection_result_list, columns=['Gene', 'Category1','Category2','Category3','Start','End','Strand','LogFC(Gene)','LogFC(5_Prime)','LogFC(3_Prime)','Qval(Gene)','Qval(5_Prime)','Qval(3_Prime)','Log_CPM(Gene)','Log_CPM(5_Prime)','Log_CPM(3_Prime)','confidence_score_upregulated','confidence_score_downregulated'])
        print("Result_DF after knockout-protection",self.result_df.head())
        # self.result_df.columns = ['Gene', 'Category1','Category2','Category3','Start','End','Strand','LogFC(Gene)','LogFC(5_Prime)','LogFC(3_Prime)','Qval(Gene)','Qval(5_Prime)','Qval(3_Prime)','Log_CPM(Gene)','Log_CPM(5_Prime)','Log_CPM(3_Prime)','confidence_score_upregulated','confidence_score_downregulated']
        self.result_df = self.result_df.dropna(how="all")
        self.result_df.to_csv("GeneReport_Post_categorize_knockout_protection.csv",index=False)
        # print("Time to category1 categorization new",time.time()-time_category1_categorization)
        misc1_time=time.time()


        # filter genes based on logfc significance (qval<5%) for forward & reverse compare insertions
        self.forward_compare_csv= self.read_and_filter_logfc_qval(self.forward_compare_csv,"prime_ends")
        self.reverse_compare_csv= self.read_and_filter_logfc_qval(self.reverse_compare_csv,"prime_ends")
        self.merged_forward_reverse_compare_csv = pd.merge(self.forward_compare_csv, self.reverse_compare_csv, on='gene_name', how='outer', suffixes=('_forward', '_reverse'))
        # print("Miscellanous activity1 filter forward-reverse, merge etc",time.time()-misc1_time)
        get_unique_gene_time= time.time()
        unique_gene_names = self.get_unique_genes_that_have_significant_3_5_prime_logfc()
        # print("Get unique list time",time.time()-get_unique_gene_time)
        up_down_categorization_time= time.time()
        for gene_name in unique_gene_names:
            #check upregulation and downregualtion of a gene based on 3,5 prime insertions
            self.prime_end_expression_categorization(gene_name)
        self.result_df.to_csv("GeneReport_Post_prime_end_expression_categorization.csv",index=False)
        # print("Time to up-downregulation categorization 1",time.time()-up_down_categorization_time)


        time_ins_count_check=time.time()
        self.result_df= self.insertion_threshold_check(self.result_df)
        # print("Time to Insertion count filtering",time.time()-time_ins_count_check)
        self.result_df.to_csv("GeneReport_Post_insertion_threshold_check.csv",index=False)

        up_down_regulation_gene_names= self.get_genes_with_up_or_down_regulation()
        self.result_df["Category4"]=None
        cat_exp_time=time.time()
        for gene in up_down_regulation_gene_names:
            #check acivity on or because of ajacent gene for upregulated or downregulate gene 
            self.categorize_activitation_explanation(gene)
        self.result_df.to_csv("GeneReport_Post_categorize_activitation_explanation.csv",index=False)
        # print("Time to categorization explanation",time.time()-cat_exp_time)
        # Replace NaN, None, and blanks
        self.result_df.replace(r'^\s*$', np.nan, regex=True, inplace=True)  # Replace blanks with NaN
        
        # final_result_df= self.insertion_threshold_check(self.result_df)
        # final_result_df.to_csv("GeneReport_Post_insertion_threshold_check.csv",index=False)
        
        cscore_time=time.time()
        self.result_df= self.signal_strength_score(self.result_df)
        # print("Confidence score check",time.time()-cscore_time)
        self.result_df.to_csv("GeneReport_Post_signal_strength_score.csv",index=False)
        # Drop rows where Category1, Category2, Category3, and Category4 are all NaN
        # final_result_df.dropna(subset=['Category1', 'Category2', 'Category3', 'Category4'], how='all', inplace=True)
        self.result_df.fillna("None", inplace=True)
        cols = ['Category1', 'Category2', 'Category3', 'Category4']
        self.result_df = self.result_df.loc[~self.result_df[cols].apply(
            lambda row: all(
                (x is None) or 
                (pd.isna(x)) or 
                (str(x).strip().lower() == 'none') 
                for x in row
            ),
            axis=1
        ),:]
        print("self.result_df",self.result_df)

        return self.result_df[['Gene', 'Category1','Category2','Category3','Category4','Start','End','Strand',
                            'LogFC(Gene)','LogFC(5_Prime)','LogFC(3_Prime)',
                            'Qval(Gene)','Qval(5_Prime)','Qval(3_Prime)',
                            'Log_CPM(Gene)','Log_CPM(5_Prime)','Log_CPM(3_Prime)',
                            'confidence_score_upregulated','confidence_score_downregulated']]
    

    def annotate_genes(self):
        genes = []
        for gene_number, f in enumerate(self.features):
            overlapping_blocks = self.blocks_overlapping_feature(f)

            if len(overlapping_blocks) == 0:
                # no hits to any blocks so move to next feature
                continue

            g = Gene(f, overlapping_blocks)

            # only consider block at a time
            for b in overlapping_blocks:
                g.upstream.append(self.find_upstream_gene(b, gene_number))
                if self.is_feature_contained_within_block(b, f):
                    g.categories.append('knockout')
                elif self.is_block_near_end_of_feature(b, f):
                    if b.max_logfc > 0.0:
                        g.categories.append('increased_mutants_at_end_of_gene')
                    else:
                        g.categories.append('decreased_mutants_at_end_of_gene')
                elif self.is_block_near_start_of_feature(b, f):
                    if b.max_logfc > 0.0:
                        g.categories.append('increased_mutants_at_start_of_gene')
                    else:
                        g.categories.append('decreased_mutants_at_start_of_gene')

            if len(g.categories) == 0:
                p = self.proportion_blocks_overlap_with_gene(f, overlapping_blocks)
                if p > 0.9:
                    g.categories.append('knockout')
                elif p > 0.8:
                    g.categories.append('knockout')
                elif p > 0.7:
                    g.categories.append('knockout')
                elif p > 0.6:
                    g.categories.append('knockout')
                elif p > 0.5:
                    g.categories.append('over_50_perc_inactivation')

            if len(g.categories) == 0:
                g.categories.append('unclassified')

            g.max_logfc_from_category()
            genes.append(g)

        # intergenic test
        intergenic_blocks = [block for block in self.blocks if block.num_genes == 0]
        for block in intergenic_blocks:
            block.upstream = self.find_nearest_upstream_gene(block)
            block.intergenic = True

        reannotate_with_5_3_prime = self.reannotate_5_3_prime(genes)

        return reannotate_with_5_3_prime

    def feature_to_gene_name(self, feature):
        gene_name_val = str(feature.location.start) + "_" + str(feature.location.end)
        if "gene" in feature.qualifiers:
            gene_name_val = feature.qualifiers["gene"][0]
        return gene_name_val

    def find_nearest_upstream_gene(self, block):
        if block.direction == 'forward':
            for f in self.features:
                if f.location.start - 1 > block.end and f.location.strand == 1:
                    return self.feature_to_gene_name(f)
        elif block.direction == 'reverse':
            for f in reversed(self.features):
                if f.location.end < block.start and f.location.strand == -1:
                    return self.feature_to_gene_name(f)
        return "NA"
    
    def feature_to_gene_name(self, feature):
        gene_name_val = str(feature.location.start) + "_" + str(feature.location.end)
        if "gene" in feature.qualifiers:
            gene_name_val = feature.qualifiers["gene"][0]
        return gene_name_val


    def find_upstream_gene(self, block, gene_number):
        if block.direction == 'reverse':
            for upstream_i in range(gene_number - 1, 0, -1):
                if self.features[upstream_i].location.strand == -1 and block.start - 1 > self.features[upstream_i].location.end:
                    return self.feature_to_gene_name(self.features[upstream_i])
        elif block.direction == 'forward':
            for upstream_i in range(gene_number + 1, len(self.features)):
                if self.features[upstream_i].location.strand == 1 and block.end < self.features[upstream_i].location.start:
                    return self.feature_to_gene_name(self.features[upstream_i])
        return "NA"

    def proportion_blocks_overlap_with_gene(self, gene, blocks):
        base_coverage = 0
        for b in blocks:
            for b_index in range(b.start - 1, b.end):
                if gene.location.start <= b_index < gene.location.end:
                    base_coverage += 1

        gene_length = gene.location.end - gene.location.start
        return base_coverage / gene_length

    def blocks_overlapping_feature(self, feature):
        overlapping_blocks = []

        for block in self.blocks:

            if (block.start - 1) > feature.location.end or feature.location.start > block.end:
                continue

            # genes are big so you are bound to hit one. Smallest in ecoli is 45bp so half it.
            for i in range(block.start - 1, block.end, 22):
                if i in feature:
                    overlapping_blocks.append(block)
                    block.num_genes += 1
                    break
        return overlapping_blocks

    def is_feature_contained_within_block(self, block, feature):
        if feature.location.start >= block.start - 1 and feature.location.end <= block.end:
            return True
        return False

    def is_block_contained_within_feature(self, block, feature):
        if block.start - 1 >= feature.location.start and block.end <= feature.location.end:
            return True
        return False

    def is_block_overlapping_feature_on_right(self, block, feature):
        if block.start - 1 < feature.location.end and block.start - 1 > feature.location.start and block.end > feature.location.end:
            return True
        return False

    def is_block_overlapping_feature_on_left(self, block, feature):
        if block.end < feature.location.end and block.end > feature.location.start and block.start - 1 < feature.location.start:
            return True
        return False
    
    def is_block_near_start_of_feature(self, block, feature):
        # forward
        if feature.location.strand == 1 and block.direction in ['forward', 'nodirection']:
            expression_end = feature.location.start + int(self.increased_expression_proportion_end * len(feature))
            if block.end <= expression_end and block.end > feature.location.start:
                return True
        # reverse
        if feature.location.strand == -1 and block.direction in ['reverse', 'nodirection']:
            expression_end = feature.location.start + int(self.increased_expression_proportion_end * len(feature))
            if block.end <= expression_end and block.end > feature.location.start:
                return True

        return False
    
    def is_block_near_end_of_feature(self, block, feature):
        # forward
        if feature.location.strand == 1 and block.direction in ['reverse', 'nodirection']:
            knock_out_start = feature.location.start + int(self.knockout_proportion_start * len(feature))
            if block.start - 1 >= knock_out_start and block.start - 1 < feature.location.end:
                return True
            
    def reannotate_5_3_prime(self, genes):
        name_to_genes = {g.gene_name: g for g in genes}

        filtered_names_to_genes = {}

        for name, gene in name_to_genes.items():
            directions = list(set([b.direction for b in gene.blocks]))

            if 'nodirection' in directions:
                continue
            # 5 prime
            res = re.search("^(.+)__([35])prime$", name)
            if res:
                found_gene_name = res.group(1)
                prime_end = res.group(2)
                if found_gene_name not in name_to_genes:
                    filtered_names_to_genes[found_gene_name] = Gene(self.embl_reader.genes_to_features[found_gene_name],
                                                                    [])
                    filtered_names_to_genes[found_gene_name].blocks = gene.blocks

                    regulation_category = self.regulation(filtered_names_to_genes[found_gene_name].feature.location.strand,
                                                          prime_end, directions)
                    if regulation_category:
                        filtered_names_to_genes[found_gene_name].categories.append(regulation_category)
                    else:
                        del filtered_names_to_genes[found_gene_name]

                else:
                    filtered_names_to_genes[found_gene_name] = name_to_genes[found_gene_name]
                    regulation_category = self.regulation(filtered_names_to_genes[found_gene_name].feature.location.strand,
                                                          prime_end, directions)
                    if regulation_category:
                        filtered_names_to_genes[found_gene_name].categories.append(regulation_category)

        # carry over non prime genes
        for name, gene in name_to_genes.items():
            res = re.search("^(.+)__[35]prime$", name)
            if not res:
                if name not in filtered_names_to_genes:
                    filtered_names_to_genes[name] = gene

        return [g for g in filtered_names_to_genes.values()]
    
    def regulation(self, strand, prime, directions):
        if prime == '5' and strand == 1 and 'forward' in directions:
            return 'upregulated'
        elif prime == '5' and strand == -1 and 'reverse' in directions:
            return 'upregulated'
        elif prime == '3' and strand == 1 and 'reverse' in directions:
            return 'downregulated'
        elif prime == '3' and strand == -1 and 'forward' in directions:
            return 'downregulated'
        else:
            return None

