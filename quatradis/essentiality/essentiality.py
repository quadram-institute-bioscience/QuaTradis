import csv
import os
import shutil
import subprocess
from pathlib import Path
from tempfile import mkstemp

import quatradis.tisp.analyse as tisp_analyse
import quatradis.util.file_handle_helpers as fhh

from quatradis.essentiality.prepare_plot_files import PrepareInputFiles

class GeneEssentiality:
    def __init__(self):
        self.condition = 0
        self.control = 0
        self.number_of_reps = 1

    def status(self):

        if self.condition == self.number_of_reps and self.control == 0:
            return 'conditionally essential'
        elif self.condition / self.number_of_reps < 0.5 and self.condition / self.number_of_reps > 0.5:
            return 'probably conditionally essential'
        elif self.condition == 0 and self.control == self.number_of_reps:
            return 'essential in control'
        elif self.control / self.number_of_reps > 0.5 and self.condition / self.number_of_reps < 0.5:
            return 'probably essential in control'
        elif self.condition == self.control and self.condition == self.number_of_reps:
            return 'always essential'
        elif (self.control / self.number_of_reps > 0.5) and (self.condition / self.number_of_reps > 0.5):
            return 'probably always essential'
        elif self.condition == 0 and self.control == 0:
            return 'always non-essential'
        elif (self.condition / self.number_of_reps < 0.5) and (self.control / self.number_of_reps < 0.5):
            return 'probably always non-essential'
        else:
            return 'inconsistent replicates'

class PlotEssentiality:
    def __init__(self, plotfile_obj, gene_insert_sites_filename, tradis_essentiality_filename, type,
                 only_essential_filename):
        self.plotfile_obj = plotfile_obj
        self.gene_insert_sites_filename = gene_insert_sites_filename
        self.tradis_essentiality_filename = tradis_essentiality_filename
        self.only_essential_filename = only_essential_filename
        self.type = type


class PlotAllEssentiality:
    def __init__(self, forward, reverse, combined, original, embl_filename):
        self.forward = forward
        self.reverse = reverse
        self.combined = combined
        self.original = original
        self.embl_filename = embl_filename





def analyse_insert_sites(embl_filename, plotfile, filetype, output_filename=None):
    if not output_filename:
        fd, output_filename = mkstemp()
    analysis_file = tisp_analyse.count_insert_sites(embl_filename, [plotfile], output_dir=os.path.dirname(output_filename), output_suffix=filetype + ".stats")
    if analysis_file != output_filename:
        shutil.move(analysis_file, output_filename)
    return output_filename


def essentiality(embl_filename, plotfile, filetype, plotname, prefix="essentiality", verbose=False):

    analysis_file = analyse_insert_sites(embl_filename, plotfile, filetype, prefix+".stats")

    e = TradisEssentialityRunner(analysis_file, verbose, prefix=prefix, analysis_type=filetype)
    e.run(plotname)
    pe = PlotEssentiality(plotfile, analysis_file, e.output_filename, filetype, e.essential_filename)

    if verbose:
        print("Essentiality:\t" + filetype + "\t" + e.output_filename)

    return pe


def run_essentiality(prepared_embl, plotfile, output_dir, minimum_threshold, plotname="", original_embl=None, verbose=False):

    if not original_embl:
        original_embl = prepared_embl

    if not plotname:
        plotname = os.path.basename(plotfile).split(sep=".")[0]

    fhh.ensure_output_dir_exists(output_dir)

    p = PrepareInputFiles(plotfile, minimum_threshold, True)

    p.create_all_files(output_dir)

    f = essentiality(prepared_embl, p.forward_plot_filename, 'forward', plotname, prefix=p.forward_plot_filename, verbose=verbose)
    r = essentiality(prepared_embl, p.reverse_plot_filename, 'reverse', plotname, prefix=p.reverse_plot_filename, verbose=verbose)
    c = essentiality(prepared_embl, p.combined_plot_filename, 'combined', plotname, prefix=p.combined_plot_filename, verbose=verbose)
    o = essentiality(original_embl, plotfile, 'original', plotname, prefix=os.path.join(output_dir, "original.plot.gz"), verbose=verbose)

    return PlotAllEssentiality(f, r, c, o, prepared_embl)


class TradisEssentialityRunner:
    def __init__(self, tabfile, verbose, exec="tradis_essentiality.R", prefix="", analysis_type=""):
        self.tabfile = tabfile
        self.exec = exec
        self.verbose = verbose
        self.prefix = prefix
        self.analysis_type = analysis_type

        self.check_if_runnable(exec)

        self.output_filename = tabfile + ".all.tsv"
        self.essential_filename = tabfile + ".essen.csv"

    def construct_command(self):
        return " ".join([self.exec, self.tabfile])

    def check_if_runnable(self, exec):

        r_script = os.path.join(os.path.dirname(__file__), exec)

        if os.path.exists(r_script):
            self.exec = r_script
        else:
            # Otherwise try to find it on the PATH
            if not shutil.which(self.exec):
                raise Exception("Can't find " + self.exec + " available to execute")

        if not shutil.which("Rscript"):
            raise Exception("Can't find Rscript executable.  Make sure R is installed and properly configured on your system and try again.")

    def run(self, plotname):
        cmd = self.construct_command()
        if self.verbose:
            print(cmd)
        subprocess.check_output(cmd, shell=True)

        self.replace_comma_tabs(self.tabfile + ".all.csv", self.output_filename)

        if self.verbose:
            print("all.csv\t" + self.output_filename)
            print("essen.csv\t" + self.essential_filename)

        if self.analysis_type == "original":
            condition_name = os.path.join(self.prefix, plotname + "." + self.analysis_type + ".ess")
            os.makedirs(Path(condition_name).parent.absolute(), exist_ok=True)
            shutil.copy(self.tabfile + ".essen.csv", condition_name)

        return self

    def replace_comma_tabs(self, input_file, output_file):
        with open(input_file, newline='') as csvfile:
            comparison_reader = csv.reader(csvfile, delimiter=',')

            with open(output_file, 'w') as outputfh:
                for line in comparison_reader:
                    outputfh.write("\t".join(line) + "\n")
        os.remove(input_file)
