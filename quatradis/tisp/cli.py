
import quatradis.tisp.generator.cli as generators

from quatradis.util.parser import create_parser

from quatradis.tisp.combine import combine
from quatradis.tisp.analyse import count_insert_sites
from quatradis.tisp.normalise import NormalisePlots

def add_subparser(subparsers):
    plot_parser_desc = "TraDIS plot file tools"
    plot_parser = subparsers.add_parser("plot", help=plot_parser_desc)
    plot_subparsers = plot_parser.add_subparsers(title=plot_parser_desc)

    generators.add_subparser(plot_subparsers)

    create_parser("combine", plot_subparsers, combine_plots, combine_plot_options,
                  "Combine multiple TraDIS plot files into one.",
                  description='''Combine multiple TraDIS plot files into one.  
    Also can combine plots files in groups.  The plot combining process is controlled via a tab-delimited file with an ID as 
    the first column followed by a list of plotfiles to combine per row. The ID will be used to name the new plotfile and as 
    an identifier in the stats file, so ensure these are unique.''')
    create_parser("count", plot_subparsers, count_plot, count_plot_options,
                  "Produce summary of insert site details from one or more plot files",
                  description="Take in one or more plot files and an embl file to produce a tab delimited file with insert site details to use as input to another script to test for essentiality.")
    create_parser("normalise", plot_subparsers, normalise_plots, normalise_plot_options,
                  "Given a set of plot files, find the one with the highest number of reads and normalise all the rest into new files",
                  usage="tradis tisp normalise [options] <plot_file>+")



def combine_plot_options(parser):
    parser.add_argument('plot_file_list', type=str,
                        help='A tab delimited file with id in the first column followed by a list of plot files')
    parser.add_argument('--combined_dir', dest='combined_dir', default="combined",
                        help='The directory in which to store the combined output files (default: combined)')


def combine_plots(args):
    combine(args.plot_file_list, args.combined_dir)


def count_plot_options(parser):
    parser.add_argument('embl_in', type=str,
                        help='The embl formatted annotation file')
    parser.add_argument('plot_in', type=str, nargs='+',
                        help='The insertion site plot files (can be gzipped)')
    parser.add_argument('-o', '--output_dir', type=str, default="",
                        help='The directory in which to put all output files (default: current working directory)')
    parser.add_argument('-s', '--output_suffix', type=str, default="tradis_gene_insert_sites.csv",
                        help='The suffix to add to output files (optional, default = tradis_gene_insert_sites.csv)')
    parser.add_argument('--trim5', action='store_true',
                        help="Trim insertion sites from 5' end of gene")
    parser.add_argument('--trim3', action='store_true',
                        help="Trim insertion sites from 3' end of gene")
    parser.add_argument('-j', '--joined_output', action='store_true',
                        help="output a single file with all info.")


def count_plot(args):
    count_insert_sites(args.embl_in, args.plot_in, args.joined_output,
                       output_dir=args.output_dir, output_suffix=args.output_suffix,
                       trim5=args.trim5, trim3=args.trim3)


def normalise_plot_options(parser):
    parser.add_argument('plot_files', type=str, nargs='+',
                        help='One or more paths to transposon insertion site plot files')
    parser.add_argument('--minimum_proportion_insertions', '-d',
                        help='If the proportion of insertions is too low compared to control, dont call decreased insertions below this level',
                        type=float, default=0.1)


def normalise_plots(args):
    n = NormalisePlots(args.plot_files, args.minimum_proportion_insertions, output_temp_files=False)
    plotfiles, max_plot_reads = n.create_normalised_files()
    n.decreased_insertion_reporting(max_plot_reads)