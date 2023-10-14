import csv
import os
import shutil
import subprocess
from pathlib import Path
from tempfile import mkstemp

import quatradis.tisp.analyse as tisp_analyse
import quatradis.util.file_handle_helpers as fhh


class GeneEssentiality:
    def __init__(self):
        self.condition = 0
        self.control = 0
        self.number_of_reps = 1

    def status(self):

        if self.condition == self.number_of_reps and self.control == 0:
            return 'conditionally essential'
        elif self.condition / self.number_of_reps < 0.5 < self.condition / self.number_of_reps:
            return 'probably conditionally essential'
        elif self.condition == 0 and self.control == self.number_of_reps:
            return 'essential in control'
        elif self.control / self.number_of_reps > 0.5 > self.condition / self.number_of_reps:
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


def essentiality(analysis_file, verbose=False):
    e = TradisEssentialityRunner(analysis_file, verbose)
    e.run()


class TradisEssentialityRunner:
    def __init__(self, count_file, verbose, exec="tradis_essentiality.R"):
        self.count_file = count_file
        self.exec = exec
        self.verbose = verbose

        self.check_if_runnable(exec)

    def construct_command(self):
        return " ".join([self.exec, self.count_file])

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

    def run(self):
        cmd = self.construct_command()
        if self.verbose:
            print(cmd)
        subprocess.check_output(cmd, shell=True)

        if self.verbose:
            print("all.csv\t" + self.count_file + ".all.csv")
            print("essen.csv\t" + self.count_file + ".essen.csv")

        return self

