import os

from quatradis.gene.block_identifier import BlockIdentifier
from quatradis.gene.annotator import GeneAnnotator
from quatradis.gene.gene import Gene

from quatradis.util.file_handle_helpers import ensure_output_dir_exists


def gene_statistics(comparison_dir, window_size, embl_file, output_dir="output", annotation_file=None):

    combined_plotfile = os.path.join(comparison_dir, "combined.logfc.plot")
    forward_plotfile = os.path.join(comparison_dir, "forward.logfc.plot")
    reverse_plotfile = os.path.join(comparison_dir, "reverse.logfc.plot")

    combined_scorefile = os.path.join(comparison_dir, "combined.pqvals.plot")

    ensure_output_dir_exists(output_dir)

    use_annotation = True if annotation_file else False

    b = BlockIdentifier(combined_plotfile, forward_plotfile, reverse_plotfile, combined_scorefile, window_size)
    blocks = b.block_generator()
    ant_file = embl_file
    if use_annotation:
        ant_file = annotation_file

    genes = GeneAnnotator(ant_file, blocks).annotate_genes()
    intergenic_blocks = [block for block in blocks if block.intergenic]

    if not use_annotation:
        all_genes = merge_windows(genes)
    else:
        all_genes = []
        for g in genes:
            all_genes.append(g)

    if len(all_genes) == 0 and len(intergenic_blocks) == 0:
        print("No significant genes found for chosen parameters.\n")
        return []

    write_gene_report(all_genes, intergenic_blocks, output_dir, use_annotation)
    write_regulated_gene_report(all_genes, output_dir)

    # if self.verbose:
    # self.print_genes_intergenic(genes,intergenic_blocks)

    return genes


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



def write_gene_report(genes, intergenic_blocks, output_dir, use_annotation):
    block_filename = os.path.join(output_dir, "gene_report.csv")

    with open(block_filename, 'w') as bf:
        bf.write(str(genes[0].header()) + "\n")
        if not use_annotation:
            for i in genes:
                bf.write(i.window_string() + "\n")
        else:
            for i in genes:
                bf.write(str(i) + "\n")
        for b in intergenic_blocks:
            bf.write(str(b) + "\n")


def write_regulated_gene_report(genes, output_dir):
    regulated_genes = [g for g in genes if g.category() == 'upregulated' or g.category() == 'downregulated']
    if len(regulated_genes) > 0:
        block_filename = os.path.join(output_dir, "regulated_gene_report.csv")
        with open(block_filename, 'w') as bf:
            bf.write(str(regulated_genes[0].header()) + "\n")
            for i in regulated_genes:
                bf.write(str(i) + "\n")


def print_genes_intergenic_blocks(genes, intergenic_blocks):
    print(genes[0].header())
    for i in genes:
        print(i)
    for b in intergenic_blocks:
        print(b)
