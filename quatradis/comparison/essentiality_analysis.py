import csv
import os
from dataclasses import dataclass

from quatradis.comparison.essentiality import GeneEssentiality
from quatradis.util.file_handle_helpers import ensure_output_dir_exists


@dataclass
class EssentialityInput:
    control_files: list
    condition_files: list
    only_ess_control_files: list
    only_ess_condition_files: list


def gene_names_from_essentiality_file(filename):
    with open(filename, "r") as fileh:
        reader = csv.reader(fileh, delimiter=",", quotechar='"')
        gene_names = [r[1] for r in reader if r[1] != "gene_name"]
        print("Number of all genes:" + str(len(gene_names)))

    return gene_names


def get_all_gene_names(control_files, condition_files):
    all_gene_names = set()

    for filename in condition_files:
        with open(filename, "r") as fileh:
            reader = csv.reader(fileh, delimiter="\t", quotechar='"')
            gene_names1 = [r[1] for r in reader if len(r) > 1 and r[1] != "gene_name"]
            all_gene_names = all_gene_names.union(set(gene_names1))

    for filename in control_files:
        with open(filename, "r") as fileh:
            reader = csv.reader(fileh, delimiter="\t", quotechar='"')
            gene_names2 = [r[1] for r in reader if len(r) > 1 and r[1] != "gene_name"]
            all_gene_names = all_gene_names.union(set(gene_names2))

    return list(all_gene_names)


def all_gene_essentiality(input: EssentialityInput, analysis_type, verbose=False):

    all_gene_names = get_all_gene_names(input.control_files, input.condition_files)
    if verbose:
        print("# all_gene_names: " + str(len(all_gene_names)))

    genes_ess = {g: GeneEssentiality() for g in all_gene_names}
    if analysis_type == "original":
        for f in input.only_ess_condition_files:
            ess_gene_names = gene_names_from_essentiality_file(f)
            if verbose:
                print("ess_gene_names condition: " + str(len(ess_gene_names)))
                print("genes_ess: " + str(len(genes_ess)))
            for e in genes_ess:
                if e in ess_gene_names:
                    genes_ess[e].condition += 1
                genes_ess[e].number_of_reps = len(input.only_ess_condition_files)
        for f in input.only_ess_control_files:
            ess_gene_names = gene_names_from_essentiality_file(f)
            if verbose:
                print("ess_gene_names control: " + str(len(ess_gene_names)))
                print("genes_ess: " + str(len(genes_ess)))
            for e in genes_ess:
                if e in ess_gene_names:
                    genes_ess[e].control += 1
                genes_ess[e].number_of_reps = len(input.only_ess_control_files)
    else:
        for e in genes_ess:
            genes_ess[e].control = 0
            genes_ess[e].condition = 0
            genes_ess[e].number_of_reps = len(input.only_ess_control_files)

    return genes_ess


def add_gene_essentiality_to_file(
    input_filename, output_filename, genes_ess, analysis_type
):
    """
    We can add information on gene essentiality to the comparison output,
    but this will not reflect full set of essential genes as the output does not contain all genes
    In order to prevent confusion this is "switched" off, but could be used if uncommented
    """
    with open(input_filename, "r") as inputfh:
        output_content = []

        reader = csv.reader(inputfh, delimiter=",", quotechar='"')
        input_content = [r for r in reader]
        # if analysis_type == "original":
        #     print("Number of cells: " + len(input_content))
        #     for i, cells in enumerate(input_content):
        #         if i == 0:
        #             cells.append("Essentiality")
        #         elif cells[1] in genes_ess and not ("3prime" in cells[1] or "5prime" in cells[1]):
        #             cells.append(genes_ess[cells[1]].status())
        #         else:
        #             cells.append('N/A')
        #         output_content.append(cells)
        # else:
        #     output_content = input_content

        output_content = input_content

        with open(output_filename, "w") as outputfh:
            for line in output_content:
                outputfh.write(",".join(line) + "\n")


def essentiality_analysis(
    input: EssentialityInput, output_dir, analysis_type, output_filename=""
):

    # out_csv = mkstemp()
    genes_ess = all_gene_essentiality(input, analysis_type)
    # add_gene_essentiality_to_file(out_csv, output_filename, genes_ess, analysis_type)
    # os.remove(out_csv)

    ensure_output_dir_exists(output_dir)
    ess = open(os.path.join(output_dir, "essentiality.csv"), "w+")
    ess.write("Gene,Essentiality,Control,Condition,Replicates\n")
    for e in genes_ess:
        ess.write(
            ",".join(
                [
                    e,
                    str(genes_ess[e].status()),
                    str(genes_ess[e].control),
                    str(genes_ess[e].condition),
                    str(genes_ess[e].number_of_reps),
                ]
            )
            + "\n"
        )
    ess.close()
