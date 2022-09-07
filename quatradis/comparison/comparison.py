import subprocess
from tempfile import mkstemp
import os
import csv
import shutil

from quatradis.comparison.plot_logfc import PlotLog
from quatradis.essentiality.essentiality import GeneEssentiality
from quatradis.util.file_handle_helpers import ensure_output_dir_exists


def prep_essentiality_pipeline_output_for_comparison(controls_essentiality_output_dirs, conditions_essentiality_output_dirs, analysis_type):

    controls_all = []
    conditions_all = []
    controls_only_ess = []
    conditions_only_ess = []

    for ctl_dir in controls_essentiality_output_dirs:
        for x in os.walk(ctl_dir):
            all = os.path.join(x[0], analysis_type + ".plot.gz.stats.all.tsv")
            essential = os.path.join(x[0], analysis_type + ".plot.gz.stats.essen.csv")

            if os.path.exists(all):
                controls_all.append(all)
                controls_only_ess.append(essential)

    for con_dir in conditions_essentiality_output_dirs:
        for x in os.walk(con_dir):
            all = os.path.join(x[0], analysis_type + ".plot.gz.stats.all.tsv")
            essential = os.path.join(x[0], analysis_type + ".plot.gz.stats.essen.csv")

            if os.path.exists(all):
                conditions_all.append(all)
                conditions_only_ess.append(essential)

    return controls_all, controls_only_ess, conditions_all, conditions_only_ess


def run_comparisons(controls_all, conditions_all, controls_only_ess, conditions_only_ess, annotations, prefix,
                    comparison_options, verbose):

    if len(controls_all) != len(controls_all):
        raise Exception("Controls and conditions must have same number of files")

    if len(controls_only_ess) != len(conditions_only_ess):
        raise Exception("Controls and conditions containing only essential genes must have same number of files")

    if len(controls_all) != len(controls_only_ess):
        raise Exception("Number of controls and conditions files must be the same for both all genes and files containing only essential genes")

    if verbose:
        print("Running comparison between condition and control and write plot files and gene report.\n")

    forward_plotfile, forward_scoresfile = generate_logfc_plot('forward', controls_all, conditions_all,
                                                               controls_only_ess, conditions_only_ess,
                                                               annotations, prefix, comparison_options, verbose)
    if verbose:
        print("Forward insertions have been compared.\n")
    reverse_plotfile, reverse_scoresfile = generate_logfc_plot('reverse', controls_all, conditions_all,
                                                               controls_only_ess, conditions_only_ess, annotations,
                                                               prefix, comparison_options, verbose)
    if verbose:
        print("Reverse insertions have been compared.\n")
    combined_plotfile, combined_scoresfile = generate_logfc_plot('combined', controls_all, conditions_all,
                                                                 controls_only_ess, conditions_only_ess, annotations,
                                                                 prefix, comparison_options, verbose)
    if verbose:
        print("Combined insertions have been compared.\n")
    original_plotfile, original_scoresfile = generate_logfc_plot('original', controls_all, conditions_all,
                                                                 controls_only_ess, conditions_only_ess, annotations,
                                                                 prefix, comparison_options, verbose)
    if verbose:
        print("Original insertions have been compared.\n")

    return forward_plotfile, reverse_plotfile, combined_plotfile, original_plotfile, \
           forward_scoresfile, reverse_scoresfile, combined_scoresfile, original_scoresfile


def generate_logfc_plot(analysis_type, controls_all, conditions_all, controls_only_ess, conditions_only_ess,
                        annotations, prefix, options, verbose):
    t = TradisComparisonRunner(conditions_all, controls_all, verbose, options.minimum_block,
                               conditions_only_ess,
                               controls_only_ess, analysis_type, prefix)
    t.run()
    p = PlotLog(t.output_filename, annotations, options,
                output_plot_filename=os.path.join(prefix, analysis_type + ".logfc.plot"),
                output_scores_filename=os.path.join(prefix, analysis_type + ".pqvals.plot"))
    p.construct_plot_file()

    # Move temp files to output dir
    renamed_csv_file = os.path.join(prefix, analysis_type + ".csv")
    shutil.move(t.output_filename, renamed_csv_file)

    if verbose:
        print("Comparison:\t" + renamed_csv_file)
        print("Plot log:\t" + p.plot_filename)
    return p.plot_filename, p.scores_filename





class TradisComparisonRunner:
    def __init__(self, condition_files, control_files, verbose, minimum_block, only_ess_files_condition,
                 only_ess_files_control, analysis_type, prefix, exec="tradis_comparison.R"):
        self.condition_files = condition_files
        self.control_files = control_files
        self.verbose = verbose
        self.minimum_block = minimum_block
        self.only_ess_files_condition = only_ess_files_condition
        self.only_ess_files_control = only_ess_files_control
        self.analysis_type = analysis_type
        self.prefix = prefix

        self.check_if_runnable(exec)

        fd, self.output_filename = mkstemp()
        fd, self.conditions_fofn = mkstemp()
        fd, self.controls_fofn = mkstemp()

    def check_if_runnable(self, exec):

        r_script = os.path.join(os.path.dirname(__file__), exec)

        if os.path.exists(r_script):
            self.exec = r_script
        else:
            # Otherwise try to find it on the PATH
            if not shutil.which(self.exec):
                raise Exception("Can't find " + self.exec + " available to execute")

        if not shutil.which("Rscript"):
            raise Exception(
                "Can't find Rscript executable.  Make sure R is installed and properly configured on your system and try again.")

    def gene_names_from_essentiality_file(self, filename):

        with open(filename, 'r') as fileh:
            reader = csv.reader(fileh, delimiter=',', quotechar='"')
            gene_names = [r[1] for r in reader if r[1] != 'gene_name']
            print("Number of all genes:" + str(len(gene_names)))

        return gene_names

    def get_all_gene_names(self):

        all_gene_names = []

        for filename in self.condition_files:
            with open(filename, 'r') as fileh:
                reader = csv.reader(fileh, delimiter='\t', quotechar='"')
                gene_names1 = [r[1] for r in reader if r[1] != 'gene_name']
                all_gene_names = list(set(all_gene_names) | set(gene_names1))

        for f in self.control_files:
            with open(filename, 'r') as fileh:
                reader = csv.reader(fileh, delimiter='\t', quotechar='"')
                gene_names2 = [r[1] for r in reader if r[1] != 'gene_name']
                all_gene_names = list(set(all_gene_names) | set(gene_names2))

        return all_gene_names

    def all_gene_essentiality(self, input_filename):

        all_gene_names = self.get_all_gene_names()
        print("all_gene_names: " + str(len(all_gene_names)))
        genes_ess = {g: GeneEssentiality() for g in all_gene_names}
        if self.analysis_type == "original":
            for f in self.only_ess_files_condition:
                ess_gene_names = self.gene_names_from_essentiality_file(f)
                print("ess_gene_names condition: " + str(len(ess_gene_names)))
                print("genes_ess: " + str(len(genes_ess)))
                for e in genes_ess:
                    if e in ess_gene_names:
                        genes_ess[e].condition += 1
                    genes_ess[e].number_of_reps = len(self.only_ess_files_condition)
            for f in self.only_ess_files_control:
                ess_gene_names = self.gene_names_from_essentiality_file(f)
                print("ess_gene_names control: " + str(len(ess_gene_names)))
                print("genes_ess: " + str(len(genes_ess)))
                for e in genes_ess:
                    if e in ess_gene_names:
                        genes_ess[e].control += 1
                    genes_ess[e].number_of_reps = len(self.only_ess_files_control)
        else:
            for e in genes_ess:
                genes_ess[e].control = 0
                genes_ess[e].condition = 0
                genes_ess[e].number_of_reps = len(self.only_ess_files_control)

        return genes_ess

    def create_fofn(self):
        with open(self.conditions_fofn, 'w') as fileh:
            if len(self.condition_files) == 1:
                fileh.write(self.condition_files[0] + "\n")
                fileh.write(self.condition_files[0] + "\n")
            else:
                for i in self.condition_files:
                    fileh.write(i + "\n")

        with open(self.controls_fofn, 'w') as fileh:
            if len(self.control_files) == 1:
                fileh.write(self.control_files[0] + "\n")
                fileh.write(self.control_files[0] + "\n")
            else:
                for i in self.control_files:
                    fileh.write(i + "\n")
        return self

    def add_gene_essentiality_to_file(self, input_filename, output_filename, genes_ess):

        # We can add information on gene essentiality to the comparison output,
        # but this will not reflect full set of essential genes as the output does not contain all genes
        # In order to prevent confusion this is "switched" off, but could be used if uncommented

        with open(input_filename, 'r') as inputfh:
            output_content = []

            reader = csv.reader(inputfh, delimiter=',', quotechar='"')
            input_content = [r for r in reader]
            # if self.analysis_type == "original":
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

            with open(output_filename, 'w') as outputfh:
                for line in output_content:
                    outputfh.write(",".join(line) + "\n")
        return self

    def construct_command(self):
        return " ".join(
            [self.exec, '-f', '-t', str(self.minimum_block), '--controls', self.controls_fofn, '--conditions',
             self.conditions_fofn, '-o', os.path.join(self.prefix, self.analysis_type + ".output.csv"), "-p", os.path.join(self.prefix, self.analysis_type + ".pdf")])

    def run(self):
        self.create_fofn()
        cmd = self.construct_command()
        if self.verbose:
            print(cmd)

        ensure_output_dir_exists(self.prefix)

        if self.analysis_type != "original":
            subprocess.check_output(cmd, shell=True)

        out_csv = os.path.join(self.prefix, self.analysis_type + ".output.csv")
        genes_ess = self.all_gene_essentiality(out_csv)

        ess = open(os.path.join(self.prefix, "Essentiality.txt"), "w+")
        ess.write("Gene, Essentiality, Control, Condition, Replicates\n")
        for e in genes_ess:
            ess.write(e + ", " + str(genes_ess[e].status()) + ", " + str(genes_ess[e].control) + ", " + str(
                genes_ess[e].condition) + ", " + str(genes_ess[e].number_of_reps) + "\n")
        ess.close()

        if self.analysis_type != "original":
            self.add_gene_essentiality_to_file(out_csv, self.output_filename, genes_ess)

            os.remove(out_csv)

        return self

