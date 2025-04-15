import os
import pandas as pd
import json
import re

def read_tsv_csv(file,file_type=None):
    if file_type =="csv":
        return pd.read_csv(file)
    else:
        return pd.read_csv(file, sep='\t')
    

def get_essentiality_changepoint(file_path, fallback=None):
    with open(file_path, 'r') as f:
        file_content = f.read()

        if not file_content.strip():
            raise ValueError(f"JSON file at {file_path} is empty.")

        # Replace invalid JSON values
        file_content = re.sub(r':\s*([-+]?)Inf', r': "\1Inf"', file_content)  # handle Inf, -Inf
        file_content = file_content.replace('NaN', '"NaN"')

        try:
            data = json.loads(file_content)
        except json.JSONDecodeError as e:
            raise ValueError(f"Invalid JSON in {file_path}: {str(e)}")

    value = data.get("essential_changepoint")

    try:
        if isinstance(value, str) and value.lower() in ['inf', '-inf', 'nan']:
            return fallback
        return float(value)
    except (TypeError, ValueError):
        return fallback


def add_essentiality_col(df, changepoint):
    if changepoint is None:
        df['essentiality'] = 'NA'
    else:
        df['essentiality'] = df['ins_index'].apply(
            lambda x: 'Essential' if x < changepoint else 'Non-Essential'
        )
    return df

def unique_gene_name(df):
    df['gene_name_clean'] = df['gene_name'].str.replace(r'__(5prime|3prime)$', '', regex=True)
    unique_genes = df['gene_name_clean'].unique().tolist()
    return unique_genes

def get_complete_report(gene_names, gene_report_df, condition_combined_count_df, 
                        merged_forward_reverse_count_df, combined_compare_df, 
                        merged_forward_reverse_compare_df):
    
    columns = ['Gene', 'Category1', 'Category2', 'Category3','Category4' ,'Start', 'End', 
               'Strand', 'LogFC(Gene)', 'LogFC(3_Prime)', 'LogFC(5_Prime)',
               'Qval(Gene)', 'Qval(5_Prime)', 'Qval(3_Prime)', 
               'Log_CPM(Gene)', 'Log_CPM(5_Prime)', 'Log_CPM(3_Prime)',
               'Read_Count(Gene)', 'Read_Count(3_Prime)', 'Read_Count(5_Prime)',
               'Ins_Index(Gene)', 'Ins_Index(3_Prime)', 'Ins_Index(5_Prime)',
               'confidence_score_upregulated','confidence_score_downregulated', 'Essentiality']
    
    rows = []

    for gene in gene_names:
        new_row = dict.fromkeys(columns)
        new_row['Gene'] = gene
        five_prime_name = f"{gene}__5prime"
        three_prime_name = f"{gene}__3prime"

        gene_row = condition_combined_count_df[condition_combined_count_df['gene_name'] == gene]

        if not gene_row.empty:
            gene_row = gene_row.iloc[0]
            new_row['Start'] = gene_row['start']
            new_row['End'] = gene_row['end']
            strand = gene_row['strand']
            new_row['Strand'] = strand
            new_row['Read_Count(Gene)'] = gene_row['read_count']
            new_row['Ins_Index(Gene)'] = gene_row['ins_index']
            new_row['Essentiality'] = gene_row.get('essentiality', None)

            gene_report_row = gene_report_df[gene_report_df["Gene"] == gene]
            if not gene_report_row.empty:
                gene_report_row = gene_report_row.iloc[0]
                new_row['Category1'] = gene_report_row.get('Category1')
                new_row['Category2'] = gene_report_row.get('Category2')
                new_row['Category3'] = gene_report_row.get('Category3')
                new_row['Category4'] = gene_report_row.get('Category4')
                new_row['confidence_score_upregulated'] = gene_report_row.get('confidence_score_upregulated')
                new_row['confidence_score_downregulated'] = gene_report_row.get('confidence_score_downregulated')

            # Determine 5' and 3' keys based on strand
            if strand == 1:
                logfc_5, logfc_3 = 'logFC_forward', 'logFC_reverse'
                qval_5, qval_3 = 'q.value_forward', 'q.value_reverse'
                logcpm_5, logcpm_3 = 'logCPM_forward', 'logCPM_reverse'
                rc_5, rc_3 = 'read_count_forward', 'read_count_reverse'
                ins_5, ins_3 = 'ins_index_forward', 'ins_index_reverse'
            else:
                logfc_5, logfc_3 = 'logFC_reverse', 'logFC_forward'
                qval_5, qval_3 = 'q.value_reverse', 'q.value_forward'
                logcpm_5, logcpm_3 = 'logCPM_reverse', 'logCPM_forward'
                rc_5, rc_3 = 'read_count_reverse', 'read_count_forward'
                ins_5, ins_3 = 'ins_index_reverse', 'ins_index_forward'

            # Pull values from merged_forward_reverse_count_df
            if five_prime_name in merged_forward_reverse_count_df.index:
                count_row_5prime = merged_forward_reverse_count_df.loc[five_prime_name]
                new_row['Read_Count(5_Prime)'] = count_row_5prime[rc_5]
                new_row['Ins_Index(5_Prime)'] = count_row_5prime[ins_5]
            if three_prime_name in merged_forward_reverse_count_df.index:
                count_row_3prime = merged_forward_reverse_count_df.loc[three_prime_name]
                new_row['Read_Count(3_Prime)'] = count_row_3prime[rc_3]
                new_row['Ins_Index(3_Prime)'] = count_row_3prime[ins_3]

            # Pull values from combined_compare_df
            if gene in combined_compare_df.index:
                compare_row = combined_compare_df.loc[gene]
                new_row['LogFC(Gene)'] = compare_row.get('logFC')
                new_row['Qval(Gene)'] = compare_row.get('q.value')
                new_row['Log_CPM(Gene)'] = compare_row.get('logCPM')

            # Pull values from merged_forward_reverse_compare_df
            if five_prime_name in merged_forward_reverse_compare_df.index:
                comp_row_5prime = merged_forward_reverse_compare_df.loc[five_prime_name]
                new_row['LogFC(5_Prime)'] = comp_row_5prime[logfc_5]
                new_row['Qval(5_Prime)'] = comp_row_5prime[qval_5]
                new_row['Log_CPM(5_Prime)'] = comp_row_5prime[logcpm_5]
            if three_prime_name in merged_forward_reverse_compare_df.index:
                comp_row_3prime = merged_forward_reverse_compare_df.loc[three_prime_name]
                new_row['LogFC(3_Prime)'] = comp_row_3prime[logfc_3]
                new_row['Qval(3_Prime)'] = comp_row_3prime[qval_3]
                new_row['Log_CPM(3_Prime)'] = comp_row_3prime[logcpm_3]

        rows.append(new_row)

    result_df = pd.DataFrame(rows, columns=columns)
    return result_df


def categorize_condition_based_essentiality(result_df, control_combined_count_df):
    control_ess_map = control_combined_count_df.set_index('gene_name')['essentiality'].to_dict()

    for idx, row in result_df.iterrows():
        gene = row['Gene']
        ess1 = row['Essentiality']
        ess2 = control_ess_map.get(gene, None)

        # Skip categorization if either value is "NA"
        if ess1 == "NA" or ess2 == "NA":
            continue

        if ess2 is not None:
            if ess1 == "Essential":
                if ess2 == "Non-Essential":
                    result_df.at[idx, 'Essentiality'] = "Conditionally-Essential"
                else:
                    result_df.at[idx, 'Essentiality'] = "Essential"
            elif ess1 == "Non-Essential":
                if ess2 == "Essential":
                    result_df.at[idx, 'Essentiality'] = "Conditionally-Non-Essential"
                else:
                    result_df.at[idx, 'Essentiality'] = "Non-Essential"

    return result_df





def generate_gene_report_ui(condition_combined_count_tsv,condition_forward_count_tsv,condition_reverse_count_tsv,
                            condition_combined_changepts_json,combined_compare_csv,forward_compare_csv,reverse_compare_csv,
                            gene_report,control_combined_count_tsv,control_combined_changepts_json,output):
    
    condition_combined_count_df = read_tsv_csv(condition_combined_count_tsv)
    condition_forward_count_df = read_tsv_csv(condition_forward_count_tsv)
    condition_reverse_count_df = read_tsv_csv(condition_reverse_count_tsv)
    merged_forward_reverse_count_df=pd.merge(condition_forward_count_df, condition_reverse_count_df, on='gene_name', how='inner', suffixes=('_forward', '_reverse'))
    condition_essen_changepoint=get_essentiality_changepoint(condition_combined_changepts_json)
    combined_compare_df =read_tsv_csv(combined_compare_csv,"csv")
    forward_compare_df =read_tsv_csv(forward_compare_csv,"csv")
    reverse_compare_df =read_tsv_csv(reverse_compare_csv,"csv")
    merged_forward_reverse_compare_df=pd.merge(forward_compare_df, reverse_compare_df, on='gene_name', how='inner', suffixes=('_forward', '_reverse'))
    gene_report_df= read_tsv_csv(gene_report)
    control_combined_count_df = read_tsv_csv(control_combined_count_tsv[0])
    control_essen_changepoint=get_essentiality_changepoint(control_combined_changepts_json[0])

    condition_combined_count_df = add_essentiality_col(condition_combined_count_df,condition_essen_changepoint)
    # condition_forward_count_df = add_essentiality_col(condition_forward_count_df,condition_essen_changepoint)
    # condition_reverse_count_df = add_essentiality_col(condition_reverse_count_df,condition_essen_changepoint)
    control_combined_count_df = add_essentiality_col(control_combined_count_df,control_essen_changepoint)

    gene_names= unique_gene_name(combined_compare_df.copy())
    merged_forward_reverse_count_df.set_index('gene_name', inplace=True)
    merged_forward_reverse_compare_df.set_index('gene_name', inplace=True)
    combined_compare_df.set_index('gene_name', inplace=True)

    result_df= get_complete_report(gene_names,gene_report_df,condition_combined_count_df,merged_forward_reverse_count_df,combined_compare_df,merged_forward_reverse_compare_df)
    result_df= categorize_condition_based_essentiality(result_df,control_combined_count_df)
    result_df.to_csv(output,index=False)