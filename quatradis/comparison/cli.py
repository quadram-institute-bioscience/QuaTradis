from quatradis.util.parser import create_parser

from quatradis.comparison.scatterplot import scatterplot_run
from quatradis.comparison.presence_absence import presence_absence_run
from quatradis.comparison.comparison import prep_essentiality_pipeline_output_for_comparison, run_comparisons, generate_logfc_plot
from quatradis.comparison.gene_stats import gene_statistics
from quatradis.comparison.plot_logfc import PlotLogOptions


def add_subparser(subparsers):
    compare_parser_desc = "Comparative analysis and visualisation of TraDIS experiments (albatradis)"
    compare_parser = subparsers.add_parser("compare", help=compare_parser_desc)
    compare_subparsers = compare_parser.add_subparsers(title=compare_parser_desc)
    create_parser("analyse", compare_subparsers, comparative_analysis, add_compare_options,
                  "Comparative analysis of TraDIS experiments whilst also predicting the impact of inserts on nearby genes",
                  usage="tradis compare analyse [options] --condition_dirs <condition_plotfiles>+ --control_dirs <control_plotfiles>+ <EMBL_file>")
    create_parser("logfc_plot", compare_subparsers, logfc_plot, add_logfc_options,
                  "Run logfc analysis for a specific plot file direction: forward, reverse, combined, original over all experiments.",
                  usage="tradis compare logfc_plot [options] --condition_dirs <condition_plotfiles>+ --control_dirs <control_plotfiles>+ <EMBL_file> <analysis_type>")
    create_parser("presence_absence", compare_subparsers, presence_absence, add_pa_options,
                  "Take in gene report files and produce a heatmap",
                  usage="tradis compare presence_absence [options] <EMBLfile> <gene_reports>+")
    create_parser("figures", compare_subparsers, scatterplot, add_scatterplot_options,
                  "Create graphics of controls vs conditions",
                  usage="tradis compare figures [options] --controls control1.plot control2.plot --conditions condition1.plot condition2.plot")
    create_parser("gene_report", compare_subparsers, gene_report, add_genereport_options,
                  "Generate gene report from comparison analysis data",
                  usage="tradis compare gene_report [options] <comparison_dir> <EMBL_file>")


def add_compare_options(parser):
    parser.add_argument('emblfile', help='Annotation file in EMBL format', type=str)
    parser.add_argument('--control_dirs',
                        help='Directories from essentiality processed control plot files (optionally gzipped). There must be an equal number of condition and control directories',
                        nargs='+', type=str)
    parser.add_argument('--condition_dirs',
                        help='Directories from essentiality processed condition plot files (optionally gzipped). There must be an equal number of condition and control directories',
                        nargs='+', type=str)

    parser.add_argument('--span_gaps', '-s', help='Span a gap if it is this multiple of a window size', type=int,
                        default=1)
    parser.add_argument('--minimum_block', '-b', help='Minimum number of reads which must be in 1 block in comparison',
                        type=int, default=100)
    parser.add_argument('--minimum_logfc', '-f', help='Minimum log fold change +/-', type=float, default=1)
    parser.add_argument('--minimum_logcpm', '-c', help='Minimum log counts per million +/-', type=float, default=8.0)
    parser.add_argument('--minimum_proportion_insertions', '-d',
                        help='If the proportion of insertions is too low compared to control, do not call decreased insertions below this level',
                        type=float, default=0.1)
    parser.add_argument('--prefix', '-o', help='Output directory prefix', type=str, default='output')
    parser.add_argument('--pvalue', '-p', help='Do not report anything above this p-value', type=float, default=0.05)
    parser.add_argument('--qvalue', '-q', help='Do not report anything above this q-value', type=float, default=0.05)
    parser.add_argument('--window_size', '-w', help='Window size', type=int, default=100)
    parser.add_argument('--dont_report_decreased_insertions', '-r', help="Whether or not to report decreased insertions", action='store_true', default=False)
    parser.add_argument('--genome_length', '-g', help="The length of the genome or sequence being analysed", type=int, required=True)

def comparative_analysis(args, help='', type=str):

    if len(args.control_dirs) != len(args.condition_dirs):
        raise "Must have equal number of control and condition essentiality directories to process"

    if len(args.control_dirs) == 0:
        raise "Must have at least one control and condition essentiality directory to process"

    controls_all, controls_only_ess, conditions_all, conditions_only_ess = prep_essentiality_pipeline_output_for_comparison(args.control_dirs, args.condition_dirs, 'combined')

    options = PlotLogOptions(
        genome_length=args.genome_length,
        minimum_logfc=args.minimum_logfc,
        maximum_pvalue=args.pvalue,
        maximum_qvalue=args.qvalue,
        minimum_logcpm=args.minimum_logcpm,
        window_size=args.window_size,
        span_gaps=args.span_gaps,
        report_decreased_insertions=not args.dont_report_decreased_insertions,
        minimum_block=args.minimum_block)

    run_comparisons(controls_all, conditions_all, controls_only_ess, conditions_only_ess, args.emblfile, args.prefix, options, args.verbose)

    gene_statistics(args.prefix, args.window_size, args.emblfile, args.prefix)



def add_logfc_options(parser):
    add_compare_options(parser)
    parser.add_argument("analysis_type", help='The type of plot file to analyse: forward, reverse, combined, original', type=str)

def logfc_plot(args):
    if len(args.control_dirs) != len(args.condition_dirs):
        raise "Must have equal number of control and condition essentiality directories to process"

    if len(args.control_dirs) == 0:
        raise "Must have at least one control and condition essentiality directory to process"

    controls_all, controls_only_ess, conditions_all, conditions_only_ess = prep_essentiality_pipeline_output_for_comparison(args.control_dirs, args.condition_dirs, args.analysis_type)

    options = PlotLogOptions(
        genome_length=args.genome_length,
        minimum_logfc=args.minimum_logfc,
        maximum_pvalue=args.pvalue,
        maximum_qvalue=args.qvalue,
        minimum_logcpm=args.minimum_logcpm,
        window_size=args.window_size,
        span_gaps=args.span_gaps,
        report_decreased_insertions=not args.dont_report_decreased_insertions,
        minimum_block=args.minimum_block)

    plotfile, scoresfile = generate_logfc_plot(args.analysis_type, controls_all, conditions_all,
                                               controls_only_ess, conditions_only_ess,
                                               args.emblfile, args.prefix, options, args.verbose)

    return plotfile, scoresfile


def add_pa_options(parser):
    parser.add_argument('emblfile', help='Annotation file in EMBL format', type=str)
    parser.add_argument('genereports', help='Gene report spreadsheets', nargs='+', type=str)
    parser.add_argument('--prefix', '-o', help='Output directory prefix', type=str, default='output')


def presence_absence(args):
    presence_absence_run(args)


def add_scatterplot_options(parser):
    parser.add_argument('--controls', '-c', help='control files (use 2 or more)', type=str, nargs='+')
    parser.add_argument('--conditions', '-d', help='condition files (use 2 or more)', type=str, nargs='+')
    parser.add_argument('--window_size', '-w', help='Window size', type=int, default=50)
    parser.add_argument('--prefix', '-o', help='Output filename prefix', type=str, default='scatter')

    parser.add_argument('--normalise', '-n', action='store_true', help='normalise the files', default=False)


def scatterplot(args):
    scatterplot_run(args)

def add_genereport_options(parser):
    parser.add_argument('comparison_dir', help='Output directory from tradis compare analyse', type=str)
    parser.add_argument('embl', help='Original EMBL file used for analysis', type=str)
    parser.add_argument('--window_size', '-w', help='Window size', type=int, default=50)
    parser.add_argument('--output_dir', '-o', help='Output filename prefix', type=str, default='gene_report.csv')
    parser.add_argument('--annotations', '-a', help='EMBL file used for annotations', type=str)


def gene_report(args):
    gene_statistics(args.comparison_dir, window_size=args.window_size, embl_file=args.embl, output_dir=args.output_dir, annotation_file=args.annotations)
