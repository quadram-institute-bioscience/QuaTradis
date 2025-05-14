import os
from pathlib import Path

from quatradis.gene.block_identifier import BlockIdentifier
from quatradis.gene.annotator import GeneAnnotator
from quatradis.gene.gene import Gene
import pandas as pd
from quatradis.util.file_handle_helpers import ensure_output_dir_exists

def write_regulated_gene_report(genes, output_filename):
    regulated_genes = [g for g in genes if g.category() == 'upregulated' or g.category() == 'downregulated']
    with open(output_filename, 'w') as bf:
        bf.write(str(Gene.header()) + "\n")
        for i in regulated_genes:
            bf.write(str(i) + "\n")

def write_gene_report_new(output_dir,genes_report):
    """
    Saves the gene report as a CSV file in the specified output directory.

    Args:
        output_dir (str): The directory where the gene report will be saved.
        genes_report (pandas.DataFrame): A DataFrame containing the annotated gene report.

    Returns:
        None: The function writes the DataFrame to a CSV file and does not return any value.

    Example:
        write_gene_report_new("results", annotated_genes_df)

    Notes:
        - The function assumes `genes_report` is a valid pandas DataFrame.
        - If the output directory does not exist, the caller must ensure it is created
          before invoking this function.
    """
    output_file_path = os.path.join(output_dir, "gene_report.tsv")
    genes_report.reset_index(drop=True).to_csv(output_file_path, sep='\t', index=False)
    return None

def merge_insertion_index(condition1_countfiles,condition2_countfiles):

    """
    Merges insertion index data from two sets of count files (condition1 and condition2) 
    for combined, forward, and reverse cases.

    This function processes and merges insertion index data from two conditions. It identifies
    the minimum insertion index for genes across the combined, forward, and reverse datasets
    and generates corresponding CSV files. The function is currently under development and 
    is intended for use in generating confidence scores for downstream applications.

    Args:
        condition1_countfiles (list of str): List of file paths for count files in condition1.
        condition2_countfiles (list of str): List of file paths for count files in condition2.

    Returns:
        tuple: A tuple containing:
            - combined_ins_idx_file (pd.DataFrame): Processed DataFrame for combined case.
            - forward_ins_idx_file (pd.DataFrame): Processed DataFrame for forward case.
            - reverse_ins_idx_file (pd.DataFrame): Processed DataFrame for reverse case.

    Notes:
        - This function is a work in progress and is aimed at supporting confidence score
          generation for genes based on insertion index data.
    """

    combined_ins_idx_file=None
    forward_ins_idx_file=None
    reverse_ins_idx_file= None

    for idx in range(len(condition1_countfiles)):
        file_name_condition1 = os.path.basename(condition1_countfiles[idx])
        file_name_condition2 = os.path.basename(condition2_countfiles[idx])
        condition1_df = pd.read_csv(condition1_countfiles[idx],sep='\t')
        condition2_df = pd.read_csv(condition2_countfiles[idx],sep='\t')
        if "combined" in file_name_condition1 and "combined" in file_name_condition2:
            merged_df = pd.merge(condition1_df[['gene_name', 'strand', 'ins_index']],
                     condition2_df[['gene_name', 'strand', 'ins_index']],
                     on='gene_name',
                     suffixes=('_df1', '_df2'))

            # Calculate minimum ins_index
            merged_df['min_ins_index'] = merged_df[['ins_index_df1', 'ins_index_df2']].min(axis=1)
            combined_ins_idx_file = merged_df[['gene_name', 'strand_df1', 'min_ins_index']]
            combined_ins_idx_file.rename(columns={'strand_df1': 'strand'}, inplace=True)
        elif "forward" in file_name_condition1 and "forward" in file_name_condition2:
            merged_df = pd.merge(condition1_df[['gene_name', 'strand', 'ins_index']],
                     condition2_df[['gene_name', 'strand', 'ins_index']],
                     on='gene_name',
                     suffixes=('_df1', '_df2'))

            # Calculate minimum ins_index
            merged_df['min_ins_index'] = merged_df[['ins_index_df1', 'ins_index_df2']].min(axis=1)
            forward_ins_idx_file = merged_df[['gene_name', 'strand_df1', 'min_ins_index']]
            forward_ins_idx_file.rename(columns={'strand_df1': 'strand'}, inplace=True)
        elif "reverse" in file_name_condition1 and "reverse" in file_name_condition2:
            merged_df = pd.merge(condition1_df[['gene_name', 'strand', 'ins_index']],
                     condition2_df[['gene_name', 'strand', 'ins_index']],
                     on='gene_name',
                     suffixes=('_df1', '_df2'))
            # Calculate minimum ins_index
            merged_df['min_ins_index'] = merged_df[['ins_index_df1', 'ins_index_df2']].min(axis=1)
            reverse_ins_idx_file = merged_df[['gene_name', 'strand_df1', 'min_ins_index']]
            reverse_ins_idx_file.rename(columns={'strand_df1': 'strand'}, inplace=True)
    combined_ins_idx_file.to_csv("combined_ins_idx_file.csv",index=False)
    forward_ins_idx_file.to_csv("forward_ins_idx_file.csv",index=False)
    reverse_ins_idx_file.to_csv("reverse_ins_idx_file.csv",index=False)
    return combined_ins_idx_file, forward_ins_idx_file, reverse_ins_idx_file

def gene_statistics_new(old_algorithm,plotfiles_all,forward_count_condition,reverse_count_condition,combined_count_condition,forward_count_control,reverse_count_control,combined_compare_csv, forward_compare_csv, reverse_compare_csv, embl_file, output_dir="output",gene_categorization_params_values=None):
    
    """
    Generates a gene report by annotating genes based on input data and conditions.

    Args:
        conditions_all (list of str): List of normalized plot files for all conditions 
            (e.g., combine.plot.gz files).
        combined_compare_csv (str): Path to the combined compare CSV file containing 
            logFC, p-values, and q-values.
        forward_csv_file (str): Path to the forward compare CSV file containing 
            logFC, p-values, and q-values.
        reverse_csv_file (str): Path to the reverse compare CSV file containing 
            logFC, p-values, and q-values.
        embl_file (str): Path to the prepared EMBL file used for analysis.
        output_dir (str, optional): Output directory for saving the gene report. 
            Defaults to "output".
        annotation_file (str, optional): Path to the EMBL file used for annotations. 
            Used if `use_annotation` is set to True.
        use_annotation (bool, optional): Flag indicating whether to use the annotation 
            file instead of the prepared EMBL file. Defaults to False.

    Returns:
        None: The function writes the annotated gene report to the specified output directory.
    
    Raises:
        Exception: If there are issues with file reading, writing, or missing inputs.
    
    Example:
        gene_statistics(
            ["condition1.combined.plot.gz", "condition2.combined.plot.gz"],
            "combined.compare.csv",
            "forward.compare.csv",
            "reverse.compare.csv",
            "prepared.embl",
            output_dir="output/",
            annotation_file="annotations.embl",
            use_annotation=True
        )
    """
    
    ensure_output_dir_exists(output_dir)
    ant_file = embl_file
    # print("Embl/Ant file being used",ant_file)
    genes_report = GeneAnnotator(old_algorithm,ant_file,None,plotfiles_all,forward_count_condition,reverse_count_condition,combined_count_condition,forward_count_control,reverse_count_control,combined_compare_csv,forward_compare_csv,reverse_compare_csv,**gene_categorization_params_values).annotate_genes_new()
    write_gene_report_new(output_dir,genes_report)
    return genes_report

def gene_statistics_old(old_algorithm,combined_plotfile, forward_plotfile, reverse_plotfile, combined_scorefile, window_size, embl_file, output_dir="output", annotation_file=None):
    ensure_output_dir_exists(output_dir)

    use_annotation = True if annotation_file else False

    b = BlockIdentifier(combined_plotfile, forward_plotfile, reverse_plotfile, combined_scorefile, window_size)
    blocks = b.block_generator()
    ant_file = embl_file
    # if use_annotation:
    #     ant_file = annotation_file

    genes = GeneAnnotator(old_algorithm,ant_file, blocks).annotate_genes()
    intergenic_blocks = [block for block in blocks if block.intergenic]

    if not use_annotation:
        all_genes = merge_windows(genes)
    else:
        all_genes = []
        for g in genes:
            all_genes.append(g)

    report_file = os.path.join(output_dir, "gene_report.tsv")

    if len(all_genes) == 0 and len(intergenic_blocks) == 0:
        print("No significant genes found for chosen parameters.\n")

    write_gene_report(all_genes, intergenic_blocks, report_file, use_annotation)
    write_regulated_gene_report(all_genes, os.path.join(output_dir, "regulated_gene_report.tsv"))

    # if self.verbose:
    # self.print_genes_intergenic(genes,intergenic_blocks)

    return genes

def write_gene_report(genes, intergenic_blocks, output_filename, use_annotation):

    with open(output_filename, 'w') as bf:
        bf.write(str(Gene.header()) + "\n")
        if not use_annotation:
            for i in genes:
                bf.write(i.window_string() + "\n")
        else:
            for i in genes:
                bf.write(str(i) + "\n")
        for b in intergenic_blocks:
            bf.write(str(b) + "\n")


def merge_windows(windows):
    # TODO Double check if the modifications here make sense.  Gene's no longer have a direct start and end field,
    # instead we get this from the feature sub structure.  Possibly we should just add a start and end feature to gene,
    # which is seperate from the feature but let's see how this goes.

    start_window = windows[0]
    i = 1
    merged_windows = []
    while i < len(windows):
        next_window = windows[i]

        if (next_window.categories[0] == "increased_mutants_at_end_of_gene" or \
                next_window.categories[0] == "increased_mutants_at_start_of_gene" or \
                next_window.categories[0] == "decreased_mutants_at_end_of_gene" or \
                next_window.categories[0] == "decreased_mutants_at_start_of_gene"):
            next_window.categories[0] = "knockout"
        if (start_window.categories[0] == "increased_mutants_at_end_of_gene" \
                or start_window.categories[0] == "increased_mutants_at_start_of_gene" or \
                start_window.categories[0] == "decreased_mutants_at_end_of_gene" or \
                start_window.categories[0] == "decreased_mutants_at_start_of_gene"):
            start_window.categories[0] = "knockout"
        category_equal = next_window.categories[0] == start_window.categories[0]
        expression_equal = next_window.expression_from_blocks() == start_window.expression_from_blocks()
        direction_equal = next_window.direction_from_blocks() == start_window.direction_from_blocks()
        if next_window.start <= start_window.end and category_equal and expression_equal and direction_equal:
            start_window.end = next_window.end
            start_window.max_logfc = max(start_window.max_logfc, next_window.max_logfc)
            start_window.min_pvalue = min(start_window.min_pvalue, next_window.min_pvalue)
            start_window.min_qvalue = min(start_window.min_qvalue, next_window.min_qvalue)
        else:
            copy = Gene(start_window.feature, start_window.blocks)
            copy.start = start_window.start
            copy.end = start_window.end
            copy.max_logfc = start_window.max_logfc
            copy.min_pvalue = start_window.min_pvalue
            copy.min_qvalue = start_window.min_qvalue
            copy.gene_name = str(copy.start) + "_" + str(copy.end)
            copy.categories.append(start_window.categories[0])
            merged_windows.append(copy)
            start_window = next_window

        i += 1

    merged_windows.append(start_window)
    return merged_windows


def print_genes_intergenic_blocks(genes, intergenic_blocks):
    print(genes[0].header())
    for i in genes:
        print(i)
    for b in intergenic_blocks:
        print(b)
