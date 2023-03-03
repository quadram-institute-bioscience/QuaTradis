import os

from quatradis.tisp.generator.create import run_tradis

from quatradis.util.parser import create_parser

from quatradis.tisp.combine import combine
from quatradis.tisp.analyse import count_insert_sites
from quatradis.tisp.normalise import NormalisePlots

def add_subparser(subparsers):
    plot_parser_desc = "TraDIS plot file tools"
    plot_parser = subparsers.add_parser("plot", help=plot_parser_desc)
    plot_subparsers = plot_parser.add_subparsers(title=plot_parser_desc)

    create_parser("create", plot_subparsers, combine_plots, combine_plot_options,
                  "Create a TraDIS plot file.",
                  description='''Can create a TraDIS plot file from either fastq or alignment file input''')
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


def create_plot_options(parser):
    parser.add_argument('fastq', type=str,
                        help='Either the fastq formatted reads for processing (can be gzipped)')
    parser.add_argument('reference', type=str,
                        help='The fasta formatted reference for processing.')
    parser.add_argument('-a', '--alignments', type=str, default="",
                        help="SAM or BAM file.  If provided then we will skip aligning the fastq to the reference.")
    parser.add_argument('-o', '--output_dir', dest='output_dir', default="",
                        help='The directory in which to put all output files (default: current working directory)')
    parser.add_argument('-p', '--output_prefix', dest='output_prefix', default="quatradis_out",
                        help='The filename prefix to use for all output files (default: quatradis_out)')
    parser.add_argument('-n', '--threads', type=int, default=1,
                        help='number of threads to use when mapping and sorting (default: 1)')
    parser.add_argument('-a', '--aligner', default="bwa",
                        help='mapping tool to use (bwa, smalt, minimap2, minimap2_long) (default: bwa)')
    parser.add_argument('-ni', '--no_ref_index', dest='no_ref_index', action='store_true',
                        help='If reference is alredy indexed used this flag to skip the indexing process.')
    parser.add_argument('-m', '--mapping_score', type=int, default=30,
                        help='mapping quality must be greater than X (Default: 30)')
    parser.add_argument('-t', '--tag', type=str, default="",
                        help='the tag to remove from fastq input (default: "")')
    parser.add_argument('-mm', '--mismatch', type=int, default=0,
                        help='number of mismatches allowed when matching tag (default: 0)')


def create_plot(args):
    output_prefix = os.path.join(args.output_dir, args.output_prefix)

    run_tradis(args.fastq, args.reference, output_prefix, alignments=args.alignments, tag=args.tag,
                      mapper=args.aligner, index=not args.no_ref_index, threads=args.threads, mismatch=args.mismatch,
                      mapping_score=args.mapping_score, verbose=args.verbose)


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