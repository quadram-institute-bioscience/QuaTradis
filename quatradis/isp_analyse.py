"""
Analyse insertion site plots
"""
import csv
import os.path

import numpy as np
from Bio import SeqIO

from quatradis import file_handle_helpers


def get_cds_locations(embl_file):
    cds_coordinates = []
    with open(embl_file, 'r') as inputHandle:
        recs = SeqIO.parse(inputHandle, "embl")
        for rec in recs:
            for feature in rec.features:
                if feature.type == "CDS":
                    cds_coordinates.append([feature.location.start, feature.location.end])
    return cds_coordinates


def get_insert_sites_from_plots(plot_files, joined_output):
    insert_sites = []
    for f in plot_files:
        insert_sites.append(np.array(read_in_plot_file(f)))

    np_insert_sites = np.array(insert_sites)

    if joined_output:
        # Combines all plot files together to give a joined insert site profile
        joined = np.sum(np_insert_sites, axis=0, keepdims=True)
        return joined

    return np_insert_sites


def read_in_plot_file(plot_file):
    inserts_per_base = []
    opener, mode = file_handle_helpers.reader_opener(plot_file)
    with opener(plot_file, mode) as input_handle:
        for base in input_handle.readlines():
            inserts = base.strip().split(" ")
            combined_insert_for_base = int(inserts[0]) + int(inserts[1])
            inserts_per_base.append(combined_insert_for_base)
    return inserts_per_base


def create_output_filenames(files, joined_output=False, output_dir="", output_suffix="isp_analysis"):
    if joined_output:
        return [os.path.join(output_dir, "joined_output." + output_suffix)]

    out_files = []
    for f in files:
        file_base = os.path.basename(os.path.basename(f))
        out_files.append(os.path.join(output_dir, file_base + "." + output_suffix))

    return out_files


def output_header():
    return ['locus_tag', 'gene_name', 'ncrna', 'start', 'end', 'strand', 'read_count', 'ins_index', 'gene_length',
            'ins_count', 'fcn']


# TODO Wondering if there is faster way of doing this... seems slow.
def is_gene_within_cds(cds_coordinates, gene_feature):
    for current_cds in cds_coordinates:
        if current_cds[0] <= gene_feature.location.start and current_cds[1] >= gene_feature.location.end:
            return True
    return False


def get_feature_id(feature):
    if 'locus_tag' in feature.qualifiers:
        feature_id = feature.qualifiers['locus_tag'][0]
    elif 'ID' in feature.qualifiers:
        feature_id = feature.qualifiers['ID'][0]
    elif 'systematic_id' in feature.qualifiers:
        feature_id = feature.qualifiers['systematic_id'][0]
    else:
        feature_id = "_".join([feature.id, str(feature.strand), str(feature.location.start), str(feature.location.end)])

    # Remove quotes from feature id
    feature_id.strip('\"')

    return feature_id


def get_gene_name(feature):
    if 'gene' in feature.qualifiers:
        gene_name = feature.qualifiers['gene'][0]
    else:
        gene_name = get_feature_id(feature)
    # Replace any non-word character
    gene_name = "".join([c for c in gene_name if c.isalnum()])
    return gene_name


def get_product_value(feature):
    product = ""
    if 'product' in feature.qualifiers:
        product = feature.qualifiers['product'][0]
    if 'pseudo' in feature.qualifiers:
        return "pseudogene"
    return '"{0}"'.format(product) if len(product) > 0 else ""


def trim_read(feature, trim5=0, trim3=0):
    """
    Trim insertion sites from start or end of gene
    Number of bases trimmed are -trim5 or -trim3 parameters multiplied by gene length.
    :param feature: A BioPython SeqFeature
    :param trim5: Multiple of gene length to trim from 5' end of feature
    :param trim3: Multiple of gene length to trim from 3' end of feature
    """

    if feature.strand == 1:
        read_start = feature.location.start + (trim5 * int(feature.location.end - feature.location.start + 1))
        read_end = feature.location.end - (trim3 * int(feature.location.end - feature.location.start + 1))
    else:
        read_start = feature.location.start + (trim3 * int(feature.location.end - feature.location.start + 1))
        read_end = feature.location.end - (trim5 * int(feature.location.end - feature.location.start + 1))

    return read_start, read_end


def relevant_feature(feature, cds_coordinates):
    """
    If this feature is a gene within a known CDS then this is not relevant for us as we we'll handle that case when we get to
    the CDS feature instead.  Also if this feature is not a CDS, polypeptide or gene type then we don't need to consider it
    """

    if feature.type == "gene" and is_gene_within_cds(cds_coordinates, feature):
        return False

    if not (feature.type == "CDS" or feature.type == "polypeptide" or feature.type == "gene"):
        return False

    return True


def count_inserts(insert_sites, read_start, read_end):
    count = 0
    inserts = 0
    for j in range(read_start, read_end + 1):
        count += insert_sites[j]
        if insert_sites[j] > 0:
            inserts += 1
    return count, inserts


def create_row(feature, insert_sites, trim5=False, trim3=False):
    feature_id = get_feature_id(feature)
    gene_name = get_gene_name(feature)
    product_value = get_product_value(feature)
    rna_value = 1 if 'ncRNA' in feature.qualifiers else 0

    # Optionally trim read based on user settings
    read_start, read_end = trim_read(feature, trim5, trim3)

    count, inserts = count_inserts(insert_sites, read_start, read_end)

    ins_index = inserts / (read_end - read_start)
    gene_length = feature.location.end - feature.location.start

    return [feature_id, gene_name, rna_value, feature.location.start + 1, feature.location.end, feature.strand, count,
            ins_index, gene_length, inserts, product_value]


def analyse_insert_sites(embl_file, plot_files, joined_output=False, output_dir="", output_suffix="tradis_gene_insert_sites.csv",
                         trim5=False, trim3=False):
    """
    Take in a plot file(s) and an embl file and produce a tab delimited file with insert site details to use as input to
    an another script to test for essentiality.
    """
    file_handle_helpers.ensure_output_dir_exists(output_dir)

    output_file_names = create_output_filenames(plot_files, joined_output, output_dir, output_suffix)
    insert_site_array = get_insert_sites_from_plots(plot_files, joined_output)
    cds_coordinates = get_cds_locations(embl_file)

    for i, insert_sites in enumerate(insert_site_array):
        output_filename = output_file_names[i]
        with open(output_filename, 'w') as output_fh:
            writer = csv.writer(output_fh, delimiter='\t', quotechar="'")
            writer.writerow(output_header())
            with open(embl_file, 'r') as inputHandle:
                recs = SeqIO.parse(inputHandle, "embl")
                for rec in recs:
                    for feature in rec.features:
                        if relevant_feature(feature, cds_coordinates):
                            row = create_row(feature, insert_sites, trim5, trim3)
                            writer.writerow(row)
