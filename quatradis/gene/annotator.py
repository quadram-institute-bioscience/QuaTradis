import re
import pandas as pd
import numpy as np

from quatradis.gene.gene import Gene
from quatradis.embl.reader import EMBLReader
from quatradis.tisp.parser import PlotParser


class GeneAnnotator:
    def __init__(self, annotation_file,conditions_all,combined_path, forward_path, reverse_path):
        self.annotation_file = annotation_file
        self.condition_files= conditions_all
        self.combined_csv= combined_path
        self.forward_csv= forward_path
        self.reverse_csv= reverse_path
        # self.blocks = self.sort_blocks_by_start_coord(blocks)
        self.transcription_threshold_fraction=0.88
        self.log_fc_threshold= 2
        self.blue_red_logfc_diff_threshold=0.4
        # self.knockout_proportion_start = 0.5
        # self.increased_expression_proportion_end = 0.3
        self.embl_reader = EMBLReader(self.annotation_file)
        # self.features = self.embl_reader.features
        self.features = self.embl_reader.genes_to_features
        self.result_df=None
        self.merged_forward_reverse=None
        # self.combined_ins_idx=insertion_index_files[0]
        # self.forward_ins_idx=insertion_index_files[1]
        # self.reverse_ins_idx=insertion_index_files[2]
        self.overlap_threshold=0.6
        self.insertion_lag= 1
        self.insertion_signal_similarity_threshold=0.9

    def sort_blocks_by_start_coord(self, blocks):
        sorted_blocks = sorted((b for b in blocks), key=lambda x: x.start)
        return sorted_blocks
        
    def get_inactivation_fraction(self, f, gene, condition_plots):
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
        gene_length = gene.end - gene.start

        # Slice the arrays for both forward and reverse strands
        forward_gene_inserts_array = [arr.forward[gene.start:gene.end] for arr in condition_plots]
        reverse_gene_inserts_array = [arr.reverse[gene.start:gene.end] for arr in condition_plots]

        # Combine into matrices
        forward_sliced_matrix = np.vstack(forward_gene_inserts_array)
        reverse_sliced_matrix = np.vstack(reverse_gene_inserts_array)

        # Helper function to calculate the result index for a given matrix
        def find_result_index(sliced_matrix):
            # Find the first index where all values are positive
            valid_indices = np.all(sliced_matrix > 0, axis=0)  # True where all arrays have values > 0
            if np.any(valid_indices):  # Check if any such index exists
                first_valid_index = np.where(valid_indices)[0][0]  # First index where condition matches
                return first_valid_index
            else:
                # Fallback: Find the first positive value index across all arrays
                first_positive_indices = [
                    np.where(arr > 0)[0][0] if np.any(arr > 0) else float("inf") for arr in sliced_matrix
                ]
                return min(first_positive_indices)

        # Find result indices for both forward and reverse strands
        forward_result_index = find_result_index(forward_sliced_matrix)
        reverse_result_index = find_result_index(reverse_sliced_matrix)

        # Get the minimum of the two result indices
        result_index = min(forward_result_index, reverse_result_index)
            # If the result is infinity, return -1
        if result_index == float("inf"):
            result_index = -1

        # Normalize the result by gene length
        return (result_index + 1) / gene_length


    def modified_regulation_categorization(self, geneName):
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
        After calling `modified_regulation_categorization("GeneA")`, the function will update the 
        regulation categories (up/downregulated) and log fold changes in `result_df` for "GeneA".
        """
        f = self.features[geneName]
        strand = f.strand

        def calculate_logfc_and_category(prime_3, prime_5, is_forward):
            nonlocal up_regulation_category, down_regulation_category,logfc_3_prime,logfc_5_prime, Mlogfc

            if not prime_3.empty:
                prime3_logfc = prime_3[logfc_key_3].values
                if not pd.isna(prime3_logfc) and prime3_logfc > self.log_fc_threshold:
                    down_regulation_category = "downregulated"
                    logfc_3_prime = prime3_logfc.item()
                    # Mlogfc = (prime3_logfc * prime3_ins_idx[prime3_ins_idx['gene_name']==three_prime_name]['min_ins_index'].values[0]) / prime_3[qval_key_3].values
            if not prime_5.empty:
                prime5_logfc = prime_5[logfc_key_5].values
                if not pd.isna(prime5_logfc) and prime5_logfc > self.log_fc_threshold:
                    up_regulation_category = "upregulated"
                    logfc_5_prime = prime5_logfc.item()
            
        # Determine keys based on strand
        logfc_key_3 = "logFC_reverse" if strand == 1 else "logFC_forward"
        logfc_key_5 = "logFC_forward" if strand == 1 else "logFC_reverse"
        qval_key_3 = "q.value_reverse" if strand == 1 else "q.value_forward"
        qval_key_5 = "q.value_forward" if strand == 1 else "q.value_reverse"
        # prime3_ins_idx = self.reverse_ins_idx if strand == 1 else self.forward_ins_idx
        # prime5_ins_idx = self.forward_ins_idx if strand == 1 else self.reverse_ins_idx

        # Fetch prime data
        five_prime_name = f"{geneName}__5prime"
        three_prime_name = f"{geneName}__3prime"
        prime_3 = self.merged_forward_reverse[self.merged_forward_reverse["gene_name"] == three_prime_name]
        prime_5 = self.merged_forward_reverse[self.merged_forward_reverse["gene_name"] == five_prime_name]

        up_regulation_category, down_regulation_category,logfc_3_prime,logfc_5_prime, Mlogfc = None, None, None,None,None
        calculate_logfc_and_category(prime_3, prime_5, strand == 1)

        if geneName in self.result_df['Gene'].values:
            self.result_df.loc[self.result_df["Gene"] == geneName, ['Category2','Category3', 'LogFC(3_Prime)', 'LogFC(5_Prime)']] = [
                up_regulation_category, down_regulation_category,logfc_3_prime, logfc_5_prime]
        else:
            new_row = {
                'Gene': geneName,
                'Category1':None,
                'Category2': up_regulation_category,
                'Category3':down_regulation_category,
                'Start': f.location.start,
                'End': f.location.end,
                'Strand': f.strand,
                'LogFC(Gene)': None,
                'LogFC(3_Prime)':logfc_3_prime,
                 'LogFC(5_Prime)':logfc_5_prime
            }
            self.result_df = pd.concat([self.result_df, pd.DataFrame([new_row])], ignore_index=True)
      

    def read_and_filter_logfc(self,csv_file_path):
        """
        Reads a CSV file containing gene expression data, filters the data to retain only the rows 
        where the 'q.value' column is less than 0.05, and returns the filtered dataframe.

        The function performs the following steps:
        - Loads the CSV file into a pandas dataframe.
        - Filters the dataframe to include only rows where the 'q.value' column has a value less than 0.05.

        Args:
        - csv_file_path (str): The path to the CSV file containing the gene expression data.

        Returns:
        - df (pandas.DataFrame): A filtered dataframe containing only the rows where the 'q.value' 
        column is less than 0.05.

        Note:
        - The filtering is performed to retain only the data with statistically significant results (based on a 
        q-value threshold of 0.05).
        """
        df= pd.read_csv(csv_file_path)
        df = df[df["q.value"]<0.05]
        return df
    
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
    

    def check_insertion_signal_pattern(self,signals, kmin, kmax):
        """
        Check if multiple signals share a similar pattern based on normalized cross-correlation,
        handling zero-variance signals.

        Parameters:
            signals (list of numpy.ndarray): List of 1D signals to compare.
            kmin (int): Minimum lag (shift) for cross-correlation.
            kmax (int): Maximum lag (shift) for cross-correlation.
            threshold (float): Correlation threshold to determine similarity (e.g., 0.8).

        Returns:
            bool: True if all valid signals have pairwise correlations above the threshold, False otherwise.
            dict: Pairwise correlation values for analysis.
        """
        num_signals = len(signals)
        correlations = {}  # To store pairwise correlations
        
        # Normalize each signal
        normalized_signals = []
        valid_indices = []
        for i, signal in enumerate(signals):
            signal_std = np.std(signal)
            if signal_std == 0:  # Handle zero-variance signal
                print(f"Signal {i+1} is invalid (zero variance). Skipping normalization.")
                normalized_signals.append(None)
            else:
                normalized_signals.append((signal - np.mean(signal)) / signal_std)
                valid_indices.append(i)
        
        # Compute pairwise correlations
        for i in range(num_signals):
            for j in range(i + 1, num_signals):
                key = f"Signal_{i+1}_vs_Signal_{j+1}"
                
                # Skip if either signal is invalid
                if normalized_signals[i] is None or normalized_signals[j] is None:
                    correlations[key] = None  # Mark as invalid correlation
                    continue

                max_corr = -np.inf  # Track the maximum correlation over the lag range
                for k in range(kmin, kmax + 1):
                    if k < 0:
                        corr = np.sum(normalized_signals[i][-k:] * normalized_signals[j][:len(signals[i]) + k]) / (len(signals[i]) + k)
                    elif k > 0:
                        corr = np.sum(normalized_signals[i][:len(signals[i]) - k] * normalized_signals[j][k:]) / (len(signals[i]) - k)
                    else:
                        corr = np.sum(normalized_signals[i] * normalized_signals[j]) / len(signals[i])
                    
                    max_corr = max(max_corr, corr)  # Keep track of the maximum correlation
                
                correlations[key] = max_corr  # Store maximum correlation for this pair
        
        # Check if all valid pairwise correlations meet the threshold
        all_similar = all(corr is not None and corr >= self.insertion_signal_similarity_threshold for corr in correlations.values())
        
        return all_similar, correlations


    def insertion_signal_check(self,f, check):
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

        condition_plots=[]
        for file in self.condition_files:
            condition_plots.append(PlotParser(file))
        
        mapping = {
        (1, "5prime"): "forward",
        (-1, "5prime"): "reverse",
        (1, "3prime"): "reverse",
        (-1, "3prime"): "forward",
        }
    
        # Get the appropriate attribute based on the mapping
        array_type = mapping.get((f.strand, check))
        if array_type:
            gene_inserts_array = [
                getattr(arr, array_type)[f.location.start:f.location.end] for arr in condition_plots
            ]
        else:
            raise ValueError("Invalid combination of strand and check")

        all_similar, pairwise_correlations= self.check_insertion_signal_pattern(gene_inserts_array,-self.insertion_lag,self.insertion_lag)
        return all_similar
    
    
    def categorize_activitation_explanation(self,geneName):
        f = self.features[geneName]
        strand = f.strand
        five_prime_name = f"{geneName}__5prime"
        three_prime_name = f"{geneName}__3prime"
        if strand == 1:
            logfc_key_3_same, logfc_key_3_opp = "logFC_reverse", "logFC_forward"
            logfc_key_5_same, logfc_key_5_opp = "logFC_forward", "logFC_reverse"
        else:
            logfc_key_3_same, logfc_key_3_opp = "logFC_forward", "logFC_reverse"
            logfc_key_5_same, logfc_key_5_opp = "logFC_reverse", "logFC_forward"

        #Check for 5prime/UpRegulation
        
        prime_3 = self.merged_forward_reverse[self.merged_forward_reverse["gene_name"] == three_prime_name]
        prime_5 = self.merged_forward_reverse[self.merged_forward_reverse["gene_name"] == five_prime_name]
        if not prime_5.empty:
            all_similar= self.insertion_signal_check(self.features[five_prime_name],"5prime")
            if all_similar:
                logfc_same = prime_5[logfc_key_5_same].values
                logfc_opp= prime_5[logfc_key_5_opp].values
                if not pd.isna(logfc_same) and not pd.isna(logfc_opp):
                    prime5_overlaps= self.embl_reader.find_overlaps(five_prime_name)
                    if len(prime5_overlaps)>0:
                        for gene , overlap_pct in prime5_overlaps:
                            if overlap_pct<self.overlap_threshold:
                                if abs(logfc_same.item()-logfc_opp.item()<self.blue_red_logfc_diff_threshold):
                                    self.result_df.loc[self.result_df['Gene'] == geneName, 'Category2'] = None

                            else:
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
                self.result_df.loc[self.result_df['Gene'] == geneName, 'Category2'] = None
        
        if not prime_3.empty:
            all_similar= self.insertion_signal_check(self.features[three_prime_name],"3prime")
            if all_similar:
                logfc_same = prime_3[logfc_key_3_same].values
                logfc_opp= prime_3[logfc_key_3_opp].values
                if not pd.isna(logfc_same) and not pd.isna(logfc_opp):
                    prime3_overlaps= self.embl_reader.find_overlaps(three_prime_name)
                    if len(prime3_overlaps)>0:
                        for gene , overlap_pct in prime3_overlaps:
                            if overlap_pct<self.overlap_threshold:
                                if abs(logfc_same.item()-logfc_opp.item()<self.blue_red_logfc_diff_threshold):
                                    self.result_df.loc[self.result_df['Gene'] == geneName, 'Category3'] = None
                            else:
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
                    

                                   

    def activation_categorization(self, gene_combined_csv):
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
        geneName=gene_combined_csv["gene_name"]
        res = re.search("^(.+)__([35])prime$", geneName)
        if res:
            return [None]*10
        else:
            condition_plots=[]
            for file in self.condition_files:
                condition_plots.append(PlotParser(file))
            f= self.features[geneName]
            gene = Gene(f)
            logfc_gene= gene_combined_csv["logFC"]
            qval_gene= gene_combined_csv["q.value"]
            
            if logfc_gene>self.log_fc_threshold:
                inactivation_fraction=self.get_inactivation_fraction(f,gene,condition_plots)
                if inactivation_fraction>=self.transcription_threshold_fraction:
                    gene.categories.append("unclassified")
                else:
                    gene.categories.append("knockout")
            elif logfc_gene<- self.log_fc_threshold:
                gene.categories.append("protection")
            else:
                gene.categories.append("unclassified")

            return [geneName,gene.categories[0],None,None,gene.start,gene.end,f.strand,logfc_gene,None,None]
    
    def get_unique_genes_that_have_significant_logfc(self):
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
        df['gene_name_cleaned'] = self.merged_forward_reverse['gene_name'].str.replace(r'(__)?(3prime|5prime)', '', regex=True)
        unique_gene_names = df['gene_name_cleaned'].dropna().unique()
        unique_gene_names = list(unique_gene_names)
        return unique_gene_names

    def annotate_genes(self):
        """
            Annotates genes based on provided data and generates a DataFrame with categorized results.

            This method performs the following steps:
            1. Reads and filters logFC (log Fold Change) values from combined, forward, and reverse compare csv files.
            2. Applies activation categorization to the combined CSV data.
            3. Merges forward and reverse compare csv files to identify unique gene names.
            4. Categorizes genes based on modified regulation and activation explanations.
            5. Cleans the resulting DataFrame, replaces missing values, and returns the annotated results.

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
        # filter genes based on logfc significance (qval<5%) for combined compare insertions
        self.combined_csv= self.read_and_filter_logfc(self.combined_csv)
        # perform activation catgeorization
        self.result_df= self.combined_csv.apply(self.activation_categorization, axis=1, result_type='expand')
        self.result_df.columns = ['Gene', 'Category1','Category2','Category3','Start','End','Strand','LogFC(Gene)','LogFC(3_Prime)', 'LogFC(5_Prime)']
        # remove genes that are uncategorized
        self.result_df = self.result_df.dropna(how="all")
        # filter genes based on logfc significance (qval<5%) for forward & reverse compare insertions
        self.forward_csv= self.read_and_filter_logfc(self.forward_csv)
        self.reverse_csv= self.read_and_filter_logfc(self.reverse_csv)
        self.merged_forward_reverse = pd.merge(self.forward_csv, self.reverse_csv, on='gene_name', how='outer', suffixes=('_forward', '_reverse'))
        unique_gene_names = self.get_unique_genes_that_have_significant_logfc()
        for gene_name in unique_gene_names:
            #check upregulation and downregualtion of a gene based on 3,5 prime insertions
            self.modified_regulation_categorization(gene_name)
        up_down_regulation_gene_names= self.get_genes_with_up_or_down_regulation()
        self.result_df["Category4"]=None
        for gene in up_down_regulation_gene_names:
            #check acivity on or because of ajacent gene for upregulated or downregulate gene 
            self.categorize_activitation_explanation(gene)
        # Replace NaN, None, and blanks
        self.result_df.replace(r'^\s*$', np.nan, regex=True, inplace=True)  # Replace blanks with NaN
        # Drop rows where Category1, Category2, Category3, and Category4 are all NaN
        self.result_df.dropna(subset=['Category1', 'Category2', 'Category3', 'Category4'], how='all', inplace=True)
        self.result_df.fillna("None", inplace=True)
        return self.result_df
    
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

    def blocks_overlapping_feature(self, feature,gene_name=None):
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

