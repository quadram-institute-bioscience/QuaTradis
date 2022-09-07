import os

from quatradis.util.parser import create_parser
from quatradis.tisp.generator.from_fastq import run_tradis
from quatradis.tisp.generator.from_alignments import plot


def add_subparser(subparsers):
    plot_create_parser_desc = "TraDIS plot file creation tools"
    plot_create_parser = subparsers.add_parser("create", help=plot_create_parser_desc)
    plot_create_subparsers = plot_create_parser.add_subparsers(title=plot_create_parser_desc)
    create_parser("from_fastq", plot_create_subparsers, create_from_fastq, from_fastq_options,
                 "From a single fastq, produces a insertion site plot file.",
                 usage="tradis tisp create from_fastq [options] <fastq> <reference>")
    create_parser("from_alignments", plot_create_subparsers, create_from_alignments, create_from_alignments_options,
                  "Create insertion site plot file from already mapped reads.",
                  usage="tradis tisp create from_alignments [options] <SAM/BAM/CRAM>")

def from_fastq_options(parser):
    parser.add_argument('fastq', type=str,
                        help='The fastq formatted reads for processing (can be gzipped).')
    parser.add_argument('reference', type=str,
                        help='The fasta formatted reference for processing.')
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


def create_from_fastq(args):

    output_prefix = os.path.join(args.output_dir, args.output_prefix)

    run_tradis(args.fastq, args.reference, output_prefix, tag=args.tag,
                      mapper=args.aligner, index=not args.no_ref_index, threads=args.threads, mismatch=args.mismatch,
                      mapping_score=args.mapping_score, verbose=args.verbose)


def create_from_alignments_options(parser):
    parser.add_argument('mapped_reads', type=str,
                        help='A SAM/BAM/CRAM file to check for tradis tags')
    parser.add_argument('-m', '--mapping_score', type=int, default=30,
                        help='mapping quality must be greater than X (Default: 30)')
    parser.add_argument('--outfile', default="tradis.plot",
                        help='Output file base name for plot (Default: tradis.plot)')


def create_from_alignments(args):
    plot(args.mapped_reads, plot_out_prefix=args.outfile, cutoff_score=args.mapping_score)

