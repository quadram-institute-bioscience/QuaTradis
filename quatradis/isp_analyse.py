"""
Analyse insertion site plots
"""
import csv
import os.path

import numpy as np
from Bio import SeqIO

from quatradis.file_handle_helpers import reader_opener


def get_cds_locations(embl_file):
    cds_coordinates = []
    with open(embl_file, 'r') as inputHandle:
        recs = SeqIO.parse(inputHandle, "embl")
        for rec in recs:
            for feature in rec.features:
                if feature.type == "CDS":
                    cds_coordinates.append([feature.start, feature.end])
    return cds_coordinates


def get_insert_sites_from_plots(plot_files, joined_output):
    insert_sites = []
    for f in plot_files:
        insert_sites.append(np.array(read_in_plot_file(f)))

    np_insert_sites = np.array(insert_sites)

    if joined_output:
        # Combines all plot files together to give a joined insert site profile
        joined = np.sum(np_insert_sites, axis=0)
        return joined

    return np_insert_sites


def read_in_plot_file(plot_file):
    inserts_per_base = []
    opener, mode = reader_opener(plot_file)
    with opener(plot_file, mode) as input_handle:
        for base in input_handle.readlines():
            inserts = base.strip().split(" ")
            combined_insert_for_base = int(inserts[0]) + int(inserts[1])
            inserts_per_base.append(combined_insert_for_base)
    return inserts_per_base


def create_output_filenames(files, joined_output, output_suffix):
    if joined_output:
        return "joined_output." + output_suffix

    outfiles = []
    for f in files:
        file_base = os.path.basename(os.path.basename(f))
        outfiles.append(file_base + "." + output_suffix)

    return outfiles


def output_header():
    return ['locus_tag', 'gene_name', 'ncrna', 'start', 'end', 'strand', 'read_count', 'ins_index', 'gene_length',
            'ins_count', 'fcn']


# TODO Wondering if there is faster way of doing this... seems unnecessarily slow.
def is_gene_within_cds(cds_coordinates, gene_feature):
    for current_coords in cds_coordinates:
        if current_coords[0] > gene_feature.start and current_coords[1] < gene_feature.end:
            return True
    return False


def get_feature_id(feature):
    #feature_id = int(random.randint(0, 10000))  #TODO Why do we do this, seems redundant?
    if feature.has_tag('locus_tag'):
        feature_id, _ = feature.get_tag_values('locus_tag')
    elif feature.has_tag('ID'):
        feature_id, _ = feature.get_tag_values('ID')
    elif feature.has_tag('systematic_id'):
        feature_id, _ = feature.get_tag_values('systematic_id')
    else:
        feature_id = "_".join([feature.seq_id, feature.strand, feature.start, feature.end])

    # Remove pipes from feature id
    feature_id.replace("|", "")

    return feature_id


def get_gene_name(feature):
    if feature.has_tag('gene'):
        gene_name, _ = feature.get_tag_values('gene')
    else:
        gene_name = get_feature_id(feature)
    #TODO Replace??? with nothing.  need to figure out what this means
    #gene_name =~ s/\W//g;
    return gene_name


def get_product_value(feature):
    product = ""
    if feature.has_tag('product'):
        product, _ = feature.get_tag_values('product')
    if feature.has_tag('pseudo'):
        return "pseudogene"
    return product


def get_rna_value(feature):
    return 1 if feature.has_tag('ncRNA') else 0


def analyse_insert_sites(embl_file, plot_files, joined_output, output_suffix, trim5, trim3):
    output_file_names = create_output_filenames(plot_files, joined_output, output_suffix)
    insert_site_array = get_insert_sites_from_plots(plot_files, joined_output)
    cds_coordinates = get_cds_locations(embl_file)

    for i, insert_sites in enumerate(insert_site_array):
        output_filename = output_file_names[i]
        with open(output_filename, 'w') as output_fh:
            writer = csv.writer(output_fh, delimiter='\t')
            writer.writerow(output_header())
            with open(embl_file, 'r') as inputHandle:
                recs = SeqIO.parse(inputHandle, "embl")
                for rec in recs:
                    for feature in rec.features:
                        if feature.type == "gene" and is_gene_within_cds(cds_coordinates, feature):
                            continue
                        if feature.type == "CDS" or feature.type == "polypeptide" or feature.type == "gene":
                            continue

                        feature_id = get_feature_id(feature)
                        gene_name = get_gene_name(feature)
                        product_value = get_product_value(feature)
                        rna_value = get_rna_value(feature)

                        if feature.strand == 1:
                            read_start = feature.start + (int(feature.end - feature.start + 1) if trim5 else 0)
                            read_end = feature.end - (int(feature.end - feature.start + 1) if trim3 else 0)
                        else:
                            read_start = feature.start + (int(feature.end - feature.start + 1) if trim3 else 0)
                            read_end = feature.end - (int(feature.end - feature.start + 1) if trim5 else 0)

                        count = 0
                        inserts = 0
                        for j in range(read_start, read_end + 1):
                            count += insert_sites[j]
                            if insert_sites[j] > 0:
                                inserts += 1
                        ins_index = inserts / (read_end - read_start + 1)

                        row = [feature_id, gene_name, rna_value, feature.start, feature.end, feature.strand, count,
                               ins_index, (feature.end - feature.start + 1), inserts, product_value]
                        writer.writerow(row)
