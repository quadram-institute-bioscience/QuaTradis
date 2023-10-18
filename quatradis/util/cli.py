import sys

from quatradis.artemis.project import artemis_project_run
from quatradis.embl.expand_genes import EMBLExpandGenes
from quatradis.gene.report_sets import gene_reports_run

from quatradis.util.parser import create_parser
from quatradis.util.mapper import index_reference, extract_sequence_names

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


def extract_names_utils_options(parser):
    parser.add_argument('fasta', help='Fasta file which is to be parsed', type=str)
    parser.add_argument('--output', '-o', help='Output file path to store the sequence names', type=str, default='sequence_names.txt')


def extract_names(args):
    extract_sequence_names(args.fasta, args.output)



