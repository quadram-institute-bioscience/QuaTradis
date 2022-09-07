import sys

from quatradis.artemis.project import artemis_project_run
from quatradis.embl.expand_genes import EMBLExpandGenes
from quatradis.gene.report_sets import gene_reports_run

from quatradis.util.parser import create_parser
from quatradis.util.mapper import index_reference, extract_sequence_names

from quatradis.embl.prepare import PrepareEMBLFile
from quatradis.essentiality.essentiality import run_essentiality

import quatradis.util.tags as tags

def add_subparser(subparsers):
    utils_parser_desc = "Miscellaneous utilities."
    utils_parser = subparsers.add_parser("utils", help=utils_parser_desc)
    utils_subparsers = utils_parser.add_subparsers(title=utils_parser_desc)

    tags.add_subparser(utils_subparsers)

    create_parser("index", utils_subparsers, index_ref, index_utils_options,
                  "Index a reference using specified alignment tool.",
                  usage="tradis utils index [options] <reference fasta file> <reference name>")
    create_parser("extract_names", utils_subparsers, extract_names, extract_names_utils_options,
                  "Extracts the sequence names from a fasta file",
                  usage="tradis utils extract_names [options] <reference fasta file>")
    create_parser("annotation", utils_subparsers, annotate_embl, annotation_utils_options,
                  "Take in an EMBL file and add flanking 3 prime and 5 prime annotation.",
                  usage='tradis utils annotation [options] <EMBLfile>')
    create_parser("artemis_project", utils_subparsers, artemis_project, artemis_project_utils_options,
                  "Create an artemis project file.",
                  usage="tradis utils artemis_project [options] <reference> <experiments_metadata.csv>")
    create_parser("gene_reports", utils_subparsers, gene_reports, gene_reports_utils_options,
                  "Manipulate gene_report.csv files, such as performing set operations",
                  usage="tradis utils gene_reports [options] <gene_report.csv>+")
    create_parser("prepare_embl", utils_subparsers, prepare_embl, prepare_embl_utils_options,
                  "Prepares an embl annotations file for comparative analysis from a plotfile.  If an existing embl file is supplied then genes in that file are expanded based on data from the plot file",
                  usage="tradis utils prepare_embl [options] <plot>")
    create_parser("essentiality", utils_subparsers, essentiality, essentiality_utils_options,
                  "Determines how essential each gene is based on the transposon insertion site plot",
                  usage="tradis utils essentiality [options] <plot>")



def index_utils_options(parser):
    parser.add_argument('reference', type=str,
                        help='The fasta formatted reference for processing.')
    parser.add_argument('refname', type=str,
                        help='The name for this reference.')
    parser.add_argument('-a', '--aligner', type=str, default="bwa",
                        help='read aligning tool to use (bwa, smalt, minimap2) (default: bwa)')
    parser.add_argument('-r', '--read_len', type=int, default=100,
                        help='the average read length used (required for smalt)')
    parser.add_argument('-sk', '--smalt_k', type=int,
                        help='the average read length used')
    parser.add_argument('-ss', '--smalt_s', type=int,
                        help='the average read length used')


def index_ref(args):
    cmd, exitcode = index_reference(args.reference, args.refname, args.read_len, mapper=args.aligner)
    if exitcode != 0:
        print("Error indexing reference.  Command used:", cmd, file=sys.stderr)


def annotation_utils_options(parser):
    parser.add_argument('emblfile', help='Annotation file in EMBL format', type=str)

    parser.add_argument('--feature_size', '-s', help='Feature size', type=int, default=198)
    parser.add_argument('--outputfile', '-o', help='Output file', type=str, default='output.embl')


def annotate_embl(args):
    eeg = EMBLExpandGenes(args.emblfile, args.feature_size)
    eeg.construct_file(args.outputfile)


def artemis_project_utils_options(parser):
    parser.add_argument('reference', help='reference EMBL file', type=str)
    parser.add_argument('experiments_metadata', help='experiments metadata spreadsheet', type=str)

    parser.add_argument('--control', '-c', help='control files (can use multiple times)', type=str, action='append')
    parser.add_argument('--outputfile', '-o', help='Output filename', type=str, default='project.properties')


def artemis_project(args):
    artemis_project_run(args)


def gene_reports_utils_options(parser):
    parser.add_argument('genereports', help='Gene report spreadsheets', nargs='+', type=str)
    parser.add_argument('--prefix', '-o', help='Output directory prefix', type=str, default='output')


def gene_reports(args):
    gene_reports_run(args)

def extract_names_utils_options(parser):
    parser.add_argument('fasta', help='Fasta file which is to be parsed', type=str)
    parser.add_argument('--output', '-o', help='Output file path to store the sequence names', type=str, default='sequence_names.txt')


def extract_names(args):
    extract_sequence_names(args.fasta, args.output)

def prepare_embl_utils_options(parser):
    parser.add_argument('plotfile', help='Transposon insertion site plot file to be used', type=str)
    parser.add_argument('--output', '-o', help='Output file path to store the prepared EMBL file', type=str, default='prepared.embl')
    parser.add_argument('--emblfile', '-e', help='If provided genes in this EMBL annotations file will expanded based on data in the plotfile.', type=str, default=None)
    parser.add_argument('--minimum_threshold', '-m',
                        help='Only include insert sites with this number or greater insertions', type=int, default=5)
    parser.add_argument('--prime_feature_size', '-z',
                        help='Feature size when adding 5/3 prime block when --use_annotation', type=int, default=198)
    parser.add_argument('--window_interval', '-l', help='Window interval', type=int, default=25)
    parser.add_argument('--window_size', '-w', help='Window size', type=int, default=100)


def prepare_embl(args):
    pef = PrepareEMBLFile(args.plotfile, args.minimum_threshold, args.window_size, args.window_interval,
                    args.prime_feature_size, args.emblfile)
    pef.create_file(args.output)

def essentiality_utils_options(parser):
    parser.add_argument('plotfile', help='Transposon insertion site plot file to be used', type=str)
    parser.add_argument('emblfile', help='EMBL file containing gene annotations', type=str)
    parser.add_argument('--original_emblfile', '-e',
                        help='If provided we will use this when analysing the original plot file, otherwise we used the mandatory embl file for all analyses.',
                        type=str, default=None)
    parser.add_argument('--output_dir', '-o', help='Output directory for generated files', type=str,
                        default='prepared')
    parser.add_argument('--minimum_threshold', '-m',
                        help='Only include insert sites with this number or greater insertions', type=int,
                        default=5)

def essentiality(args):
    run_essentiality(args.emblfile, args.plotfile, args.output_dir, args.minimum_threshold, original_embl=args.original_emblfile, verbose=args.verbose)
