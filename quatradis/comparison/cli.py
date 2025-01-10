from quatradis.util.parser import create_parser

from quatradis.comparison.scatterplot import scatterplot_run
from quatradis.comparison.presence_absence import presence_absence_run
from quatradis.comparison.comparison import (
    prep_essentiality_pipeline_output_for_comparison,
    run_comparisons,
    insertion_site_comparison,
)
from quatradis.comparison.gene_stats import gene_statistics
from quatradis.comparison.plot_logfc import PlotLogOptions
from quatradis.comparison.split import split_plot as sp
from quatradis.comparison.essentiality import essentiality as ess
from quatradis.embl.prepare import PrepareEMBLFile
from quatradis.comparison.essentiality_analysis import (
    essentiality_analysis as ess_an,
    EssentialityInput,
)


def add_subparser(subparsers):
    compare_parser_desc = (
        "Comparative analysis and visualisation of TraDIS experiments (albatradis)"
    )
    compare_parser = subparsers.add_parser("compare", help=compare_parser_desc)
    compare_subparsers = compare_parser.add_subparsers(title=compare_parser_desc)
    create_parser(
        "insertion_sites",
        compare_subparsers,
        insertion_site_comparison_cli,
        add_insertion_site_comparison_options,
        "Run an insertion site comparison for a specific plot file direction: forward, reverse, combined, original over all experiments to produce a file including edgeR output (logfc, p and q values)",
        usage="tradis compare insertion_sites [options] --condition_dirs <condition_plotfiles>+ --control_dirs <control_plotfiles>+ <EMBL_file> <analysis_type>",
    )
    create_parser(
        "presence_absence",
        compare_subparsers,
        presence_absence,
        add_pa_options,
        "Take in gene report files and produce a heatmap",
        usage="tradis compare presence_absence [options] <EMBLfile> <gene_reports>+",
    )
    create_parser(
        "figures",
        compare_subparsers,
        figures,
        figures_options,
        "Create graphics of controls vs conditions",
        usage="tradis compare figures [options] --controls control1.plot control2.plot --conditions condition1.plot condition2.plot",
    )
    create_parser(
        "gene_report",
        compare_subparsers,
        gene_report,
        add_genereport_options,
        "Generate gene report from comparison analysis data",
        usage="tradis compare gene_report [options] <comparison_dir> <EMBL_file>",
    )
    create_parser(
        "split",
        compare_subparsers,
        split_plot,
        split_plot_options,
        "Splits a plot file into forward, reverse and combined files",
        usage="tradis compare split [options] <EMBL_file> <plot_file>",
    )
    create_parser(
        "essentiality",
        compare_subparsers,
        essentiality,
        essentiality_utils_options,
        "Determines how essential each gene is based on the transposon insertion site plot counts",
        usage="tradis compare essentiality [options] <plot_count_file>",
    )
    create_parser(
        "prepare_embl",
        compare_subparsers,
        prepare_embl,
        prepare_embl_utils_options,
        "Prepares an embl annotations file for comparative analysis from a plotfile.  If an existing embl file is supplied then genes in that file are expanded based on data from the plot file",
        usage="tradis compare prepare_embl [options] <plot>",
    )
    create_parser(
        "essentiality_analysis",
        compare_subparsers,
        essentiality_analysis,
        essentiality_analysis_options,
        "Compares essentiality across condition and control samples",
        usage="tradis compare essentiality_analysis [options] --controls <control_gene_files> --conditions <condition_gene_files> --ess_controls <control_essential_gene_files> --ess_conditions <condition_essential_gene_files>",
    )


def add_insertion_site_comparison_options(parser):
    parser.add_argument("emblfile", help="Annotation file in EMBL format", type=str)
    parser.add_argument(
        "analysis_type",
        help="Analysis type: forward, reverse, combined, original",
        type=str,
    )
    parser.add_argument(
        "--controls",
        help="Control plot files (optionally gzipped). There must be an equal number of condition and control files.",
        nargs="+",
        type=str,
    )
    parser.add_argument(
        "--conditions",
        help="Condition plot files (optionally gzipped). There must be an equal number of condition and control files.",
        nargs="+",
        type=str,
    )

    parser.add_argument(
        "--span_gaps",
        "-s",
        help="Span a gap if it is this multiple of a window size (default: 1)",
        type=int,
        default=1,
    )
    parser.add_argument(
        "--minimum_block",
        "-b",
        help="Minimum number of reads which must be in 1 block in comparison (default: 100)",
        type=int,
        default=100,
    )
    parser.add_argument(
        "--minimum_logfc",
        "-f",
        help="Minimum log fold change +/- (default: 1)",
        type=float,
        default=1,
    )
    parser.add_argument(
        "--minimum_logcpm",
        "-c",
        help="Minimum log counts per million +/- (default: 8.0)",
        type=float,
        default=8.0,
    )
    parser.add_argument(
        "--minimum_proportion_insertions",
        "-d",
        help="If the proportion of insertions is too low compared to control, do not call decreased insertions below this level (default: 0.1)",
        type=float,
        default=0.1,
    )
    parser.add_argument(
        "--prefix", "-o", help="Output directory prefix", type=str, default="output"
    )
    parser.add_argument(
        "--pvalue",
        "-p",
        help="Do not report anything above this p-value (default: 0.05)",
        type=float,
        default=0.05,
    )
    parser.add_argument(
        "--qvalue",
        "-q",
        help="Do not report anything above this q-value (default: 0.05)",
        type=float,
        default=0.05,
    )
    parser.add_argument(
        "--window_size", "-w", help="Window size (default: 100)", type=int, default=100
    )
    parser.add_argument(
        "--dont_report_decreased_insertions",
        "-r",
        help="Whether or not to report decreased insertions",
        action="store_true",
        default=False,
    )
    parser.add_argument(
        "--genome_length",
        "-g",
        help="The length of the genome or sequence being analysed",
        type=int,
        required=True,
    )


def insertion_site_comparison_cli(args):
    if len(args.controls) != len(args.conditions):
        raise "Must have equal number of control and condition files to process"

    if len(args.controls) == 0:
        raise "Must have at least one control and condition file to process"

    # controls_all, controls_only_ess, conditions_all, conditions_only_ess = prep_essentiality_pipeline_output_for_comparison(
    #    args.control_dirs, args.condition_dirs, args.analysis_type)

    options = PlotLogOptions(
        genome_length=args.genome_length,
        minimum_logfc=args.minimum_logfc,
        maximum_pvalue=args.pvalue,
        maximum_qvalue=args.qvalue,
        minimum_logcpm=args.minimum_logcpm,
        window_size=args.window_size,
        span_gaps=args.span_gaps,
        report_decreased_insertions=not args.dont_report_decreased_insertions,
        minimum_block=args.minimum_block,
    )
    print(args.analysis_type,"||",args.controls,"||",args.conditions,"||",args.emblfile,args.prefix,options,args.verbose)

    plotfile, scoresfile = insertion_site_comparison(
        args.analysis_type,
        args.controls,
        args.conditions,
        args.emblfile,
        args.prefix,
        options,
        args.verbose,
    )

    return plotfile, scoresfile


def add_pa_options(parser):
    parser.add_argument("emblfile", help="Annotation file in EMBL format", type=str)
    parser.add_argument(
        "genereports", help="Gene report spreadsheets", nargs="+", type=str
    )
    parser.add_argument(
        "--prefix", "-o", help="Output directory prefix", type=str, default="output"
    )


def presence_absence(args):
    presence_absence_run(args)


def figures_options(parser):
    parser.add_argument(
        "--controls", "-c", help="control files (use 2 or more)", type=str, nargs="+"
    )
    parser.add_argument(
        "--conditions",
        "-d",
        help="condition files (use 2 or more)",
        type=str,
        nargs="+",
    )
    parser.add_argument(
        "--window_size", "-w", help="Window size (default: 50)", type=int, default=50
    )
    parser.add_argument(
        "--prefix", "-p", help="Output prefix", type=str, default="figures"
    )
    parser.add_argument(
        "--normalise",
        "-n",
        action="store_true",
        help="normalise the files",
        default=False,
    )


def figures(args):
    scatterplot_run(args)


def add_genereport_options(parser):
    parser.add_argument(
        "--conditionfiles",
        help="Normalized plot files for all conditions (combine.plot.gz).",
        nargs='+',  # Accept one or more arguments as a list
        type=str,
    )
    parser.add_argument("--embl", help="Prepared EMBL file used for analysis", type=str)
    parser.add_argument("--combined_compare", help="Combined compare csv with logfc, p and q values.", type=str)
    parser.add_argument("--forward_compare", help="Forward compare csv with logfc, p and q values.", type=str)
    parser.add_argument("--reverse_compare", help="Reverse compare csv with logfc, p and q values.", type=str)
    parser.add_argument(
        "--output_dir", "-o", help="Output filename prefix", type=str, default="."
    )
    parser.add_argument(
        "--annotations", "-a", help="EMBL file used for annotations", type=str
    )
    parser.add_argument("--use_annotation", "-ua", help="Flag to use annotation file in place of prepared embl file (default: False)", action='store_true')


def gene_report(args):
    """
    Wrapper function to generate a gene report by passing command-line arguments 
    to the `gene_statistics` function.

    Args:
        args (argparse.Namespace): Parsed command-line arguments containing the following:
            - conditionfiles (list of str): Normalized plot files for all conditions (e.g., combine.plot.gz).
            - combined_compare (str): Path to the combined compare CSV file with logFC, p-values, and q-values.
            - forward_compare (str): Path to the forward compare CSV file with logFC, p-values, and q-values.
            - reverse_compare (str): Path to the reverse compare CSV file with logFC, p-values, and q-values.
            - embl (str): Path to the prepared EMBL file used for analysis.
            - output_dir (str): Output filename prefix for saving the report.
            - annotations (str, optional): Path to the EMBL file used for annotations.
            - use_annotation (bool, optional): Flag to use the annotation file instead of the EMBL file.

    Returns:
        None: The function calls `gene_statistics` and does not return a value.

    Example:
        parser = argparse.ArgumentParser()
        # Add arguments to parser as described above
        args = parser.parse_args()
        gene_report(args)
    """

    gene_statistics(
    conditions_all=args.conditionfiles,
    combined_csv_file=args.combined_compare,
    forward_csv_file=args.forward_compare,
    reverse_csv_file=args.reverse_compare,
    embl_file=args.embl,
    output_dir=args.output_dir,
    annotation_file=args.annotations,
    use_annotation=args.use_annotation,
    )


def split_plot_options(parser):
    parser.add_argument("plot_file", help="Plot file to split", type=str)
    parser.add_argument(
        "--minimum_threshold",
        "-m",
        help="Minimum number of insertions required to pass filter (default: 5)",
        type=int,
        default=5,
    )
    parser.add_argument(
        "--output_dir", "-o", help="Output directory", type=str, required=True
    )
    parser.add_argument(
        "--gzipped",
        "-g",
        help="Whether to output gzipped plot files",
        action="store_true",
    )


def split_plot(args):
    sp(
        args.plot_file,
        output_dir=args.output_dir,
        minimum_threshold=args.minimum_threshold,
        gzipped=args.gzipped,
    )


def essentiality_utils_options(parser):
    parser.add_argument(
        "count_file",
        help='Counts from a transposon insertion site plot file to be used.  i.e. generated from "tradis plot count"',
        type=str,
    )


def essentiality(args):
    ess(args.count_file, verbose=args.verbose)


def prepare_embl_utils_options(parser):
    parser.add_argument(
        "--plotfile",
        help="Transposon insertion site plot files to be used (accepts multiple paths)",
        nargs='+',  # Accept one or more arguments as a list
        type=str,
    )
    parser.add_argument(
        "--output",
        "-o",
        help="Output file path to store the prepared EMBL file",
        type=str,
        default="prepared.embl",
    )
    parser.add_argument(
        "--emblfile",
        "-e",
        help="If provided genes in this EMBL annotations file will expanded based on data in the plotfile.",
        type=str,
        default=None,
    )
    parser.add_argument(
        "--minimum_threshold",
        "-m",
        help="Only include insert sites with this number or greater insertions (default: 5)",
        type=int,
        default=5,
    )
    parser.add_argument(
        "--prime_feature_size",
        "-z",
        help="Feature size when adding 5/3 prime block when --use_annotation (default: 198)",
        type=int,
        default=400,
    )
    parser.add_argument(
        "--window_interval",
        "-l",
        help="Window interval (default: 25)",
        type=int,
        default=25,
    )
    parser.add_argument(
        "--window_size", "-w", help="Window size (default: 100)", type=int, default=100
    )
    # Modification 2
    parser.add_argument("--dynamic_window", "-dw", help="Dynamic Window for 3,5 Prime_Features (default: True)", action='store_true')


def prepare_embl(args):
    # Modification 3
    pef = PrepareEMBLFile(
        args.plotfile,
        args.minimum_threshold,
        args.window_size,
        args.window_interval,
        args.prime_feature_size,
        args.emblfile,
        args.dynamic_window
    )
    pef.create_file(args.output)


def essentiality_analysis(args):

    i = EssentialityInput(
        args.controls, args.conditions, args.ess_controls, args.ess_conditions
    )
    ess_an(i, args.output_dir, args.type)


def essentiality_analysis_options(parser):
    parser.add_argument(
        "--controls", "-c", help="control files (use 2 or more)", type=str, nargs="+"
    )
    parser.add_argument(
        "--conditions",
        "-d",
        help="condition files (use 2 or more)",
        type=str,
        nargs="+",
    )
    parser.add_argument(
        "--ess_controls",
        "-e",
        help="control files (use 2 or more)",
        type=str,
        nargs="+",
    )
    parser.add_argument(
        "--ess_conditions",
        "-f",
        help="condition files (use 2 or more)",
        type=str,
        nargs="+",
    )
    parser.add_argument(
        "--output_dir", "-o", help="Output directory", type=str, default="."
    )
    parser.add_argument(
        "--type",
        "-t",
        help="Analysis type: original, combined, forward, reverse",
        type=str,
        required=True,
    )

