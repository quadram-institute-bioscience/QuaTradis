import os
import shutil
import subprocess
from tempfile import mkstemp

from quatradis.comparison.plot_logfc import PlotLog
from quatradis.util.file_handle_helpers import ensure_output_dir_exists


def prep_essentiality_pipeline_output_for_comparison(controls_essentiality_output_dirs, conditions_essentiality_output_dirs, analysis_type):

    controls_all = []
    conditions_all = []
    controls_only_ess = []
    conditions_only_ess = []

    for ctl_dir in controls_essentiality_output_dirs:
        for x in os.walk(ctl_dir):
            all = os.path.join(x[0], analysis_type + ".plot.gz.stats.all.csv")
            essential = os.path.join(x[0], analysis_type + ".plot.gz.stats.essen.csv")

            if os.path.exists(all):
                controls_all.append(all)
                controls_only_ess.append(essential)

    for con_dir in conditions_essentiality_output_dirs:
        for x in os.walk(con_dir):
            all = os.path.join(x[0], analysis_type + ".plot.gz.stats.all.csv")
            essential = os.path.join(x[0], analysis_type + ".plot.gz.stats.essen.csv")

            if os.path.exists(all):
                conditions_all.append(all)
                conditions_only_ess.append(essential)

    return controls_all, controls_only_ess, conditions_all, conditions_only_ess


def run_comparisons(controls_all, conditions_all, annotations, prefix, comparison_options, verbose):

    if len(controls_all) != len(controls_all):
        raise Exception("Controls and conditions must have same number of files")

    if verbose:
        print("Running comparison between condition and control and write plot files and gene report.\n")

    forward_plotfile, forward_scoresfile = insertion_site_comparison('forward', controls_all, conditions_all, annotations,
                                                                     prefix, comparison_options, verbose)
    if verbose:
        print("Forward insertions have been compared.\n")
    reverse_plotfile, reverse_scoresfile = insertion_site_comparison('reverse', controls_all, conditions_all, annotations,
                                                                     prefix, comparison_options, verbose)
    if verbose:
        print("Reverse insertions have been compared.\n")
    combined_plotfile, combined_scoresfile = insertion_site_comparison('combined', controls_all, conditions_all, annotations,
                                                                       prefix, comparison_options, verbose)
    if verbose:
        print("Combined insertions have been compared.\n")
    original_plotfile, original_scoresfile = insertion_site_comparison('original', controls_all, conditions_all, annotations,
                                                                       prefix, comparison_options, verbose)
    if verbose:
        print("Original insertions have been compared.\n")

    return forward_plotfile, reverse_plotfile, combined_plotfile, original_plotfile, \
           forward_scoresfile, reverse_scoresfile, combined_scoresfile, original_scoresfile


def insertion_site_comparison(analysis_type, controls_all, conditions_all,
                              annotations, prefix, options, verbose):

    # Create CSV report file containing LogFC and PValue information for genes based on their insertion sites
    t = TradisComparisonRunner(conditions_all, controls_all, verbose, options.minimum_block,
                               analysis_type, prefix)
    t.run()

    # Create plot sytle files containing the logfc and p and q values at every position in the reference
    p = PlotLog(t.get_report_file_path(), annotations, options,
                output_plot_filename=os.path.join(prefix, analysis_type + ".logfc.plot"),
                output_scores_filename=os.path.join(prefix, analysis_type + ".pqvals.plot"))
    p.construct_plot_file()

    if verbose:
        print("Comparison:\t" + t.get_report_file_path())
        print("Plot log:\t" + p.plot_filename)
    return p.plot_filename, p.scores_filename



class TradisComparisonRunner:
    def __init__(self, condition_files, control_files, verbose, minimum_block, analysis_type, prefix, exec="tradis_comparison.R"):
        self.condition_files = condition_files
        self.control_files = control_files
        self.verbose = verbose
        self.minimum_block = minimum_block
        self.analysis_type = analysis_type
        self.prefix = prefix

        self.check_if_runnable(exec)

        fd, self.conditions_fofn = mkstemp()
        fd, self.controls_fofn = mkstemp()

    def get_report_file_path(self):
        return os.path.join(self.prefix, self.analysis_type + ".compare.csv")

    def get_volcano_file_path(self):
        return os.path.join(self.prefix, self.analysis_type + ".compare.pdf")

    def check_if_runnable(self, exec):

        r_script = os.path.join(os.path.dirname(__file__), exec)

        if os.path.exists(r_script):
            self.exec = r_script
        else:
            # Otherwise try to find it on the PATH
            if not shutil.which(exec):
                raise Exception("Can't find " + exec + " available to execute")
            else:
                self.exec = exec

        if not shutil.which("Rscript"):
            raise Exception(
                "Can't find Rscript executable.  Make sure R is installed and properly configured on your system and try again.")


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

    def construct_command(self):
        return " ".join(
            [self.exec, '-f', '-t', str(self.minimum_block), '--controls', self.controls_fofn, '--conditions',
             self.conditions_fofn, '-o', self.get_report_file_path(), "-p", os.path.join(self.prefix, self.analysis_type + ".compare.pdf")])

    def run(self):
        self.create_fofn()
        cmd = self.construct_command()
        if self.verbose:
            print(cmd)

        ensure_output_dir_exists(self.prefix)

        subprocess.check_output(cmd, shell=True)

        return self

