from quatradis.util.parser import create_parser

from quatradis.comparison.scatterplot import scatterplot_run
from quatradis.comparison.presence_absence import presence_absence_run
from quatradis.comparison.comparison import (
    prep_essentiality_pipeline_output_for_comparison,
    run_comparisons,
    insertion_site_comparison,
)
from quatradis.comparison.gene_stats import gene_statistics_new, gene_statistics_old
from quatradis.comparison.plot_logfc import PlotLogOptions
from quatradis.comparison.split import split_plot as sp
from quatradis.comparison.essentiality import essentiality as ess
from quatradis.embl.prepare import PrepareEMBLFile
from quatradis.comparison.essentiality_analysis import (
    essentiality_analysis as ess_an,
    EssentialityInput,
)
from quatradis.util.config_defaults import DYNAMIC_WINDOW_PARAMS, GENE_REPORT_PARAMS , DYNAMIC_WINDOW_HELP, GENE_REPORT_HELP
from quatradis.util.ui import generate_gene_report_ui


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
        "Generates gene report from comparison analysis data.",
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
    create_parser(
        "gene_report_ui",
        compare_subparsers,
        gene_report_ui,
        gene_report_ui_options,
        "Compares essentiality across condition and control samples, combines with gene report and gives a consolidated report for visualization onto UI.",
        usage="tradis compare gene_report_ui [options] --combined_count_tsv <combined_count_tsv for condition replicate> --forward_count_tsv <forward_count_tsv for condition replicate> --reverse_count_tsv <reverse_count_tsv for condition replicate> --combined_compare <combined_compare> --forward_compare <forward_compare> --reverse_compare <reverse_compare> --gene_report <gene_report>",
    )
def gene_report_ui_options(parser):
    """
    Defines CLI arguments for the `gene_report_ui` command used in the TraDIS comparison pipeline.

    This command compares essentiality across a condition and control using insertion data 
    (combined, forward, reverse counts), changepoint data, and a gene report. It generates 
    a consolidated output report suitable for visualization in a UI.

    Arguments added to the parser:
    - --condition_combined_count: Combined count TSV of the condition replicate.
    - --condition_forward_count: Forward count TSV of the condition replicate.
    - --condition_reverse_count: Reverse count TSV of the condition replicate.
    - --condition_combined_changepts_json: Changepoint JSON of the condition replicate.
    - --combined_compare_csv: Combined compare CSV from the comparison folder.
    - --forward_compare_csv: Forward compare CSV from the comparison folder.
    - --reverse_compare_csv: Reverse compare CSV from the comparison folder.
    - --gene_report: Gene report generated from gene_stats rule.
    - --controls_combined_count_all: Combined count TSVs for all control replicates (only first is used).
    - --controls_combined_changepts_json_all: Changepoint JSONs for all control replicates (only first is used).
    - --output / -o: Output file path for the generated UI-compatible report.
    """

    parser.add_argument("--condition_combined_count", help="Combined count tsv of condition replicate.", type=str)
    parser.add_argument("--condition_forward_count", help="Forward count tsv of condition replicate.", type=str)
    parser.add_argument("--condition_reverse_count", help="Reverse count tsv of condition replicate.", type=str)
    parser.add_argument("--condition_combined_changepts_json", help="Changepoint json of condition replicate.", type=str)
    parser.add_argument("--combined_compare_csv", help="Combined compare csv in comparison folder.", type=str)
    parser.add_argument("--forward_compare_csv", help="Forward compare csv in comparison folder.", type=str)
    parser.add_argument("--reverse_compare_csv", help="Reverse compare csv in comparison folder.", type=str)
    parser.add_argument("--gene_report", help="Gene report for categorization (output of gene_stats rule).", type=str)
    parser.add_argument("--controls_combined_count_all", help="Combined count tsvs for all control replicates(Only first one will be utilized)", type=str,nargs="+")
    parser.add_argument("--controls_combined_changepts_json_all", help="Changepoint json for all one control replicates(Only first one will be utilized).", type=str,nargs="+")
    parser.add_argument("--output", "-o", help="Output file", type=str, required=True)

def gene_report_ui(args):
    """
    Executes the `gene_report_ui` comparison step using parsed CLI arguments.

    Calls `generate_gene_report_ui()` with provided condition and control datasets,
    changepoint data, comparison results, and gene report to produce a consolidated
    output for UI visualization.
    """

    generate_gene_report_ui(condition_combined_count_tsv=args.condition_combined_count,
                            condition_forward_count_tsv=args.condition_forward_count,
                            condition_reverse_count_tsv=args.condition_reverse_count,
                            condition_combined_changepts_json=args.condition_combined_changepts_json,
                            combined_compare_csv=args.combined_compare_csv,
                            forward_compare_csv=args.forward_compare_csv,
                            reverse_compare_csv=args.reverse_compare_csv,
                            gene_report=args.gene_report,
                            control_combined_count_tsv=args.controls_combined_count_all,
                            control_combined_changepts_json=args.controls_combined_changepts_json_all,
                            output=args.output)



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
    # Old Algo Params
    parser.add_argument("--combined", help="Combined plotfile", type=str)
    parser.add_argument("--forward", help="Forward plotfile", type=str)
    parser.add_argument("--reverse", help="Reverse plotfile", type=str)
    parser.add_argument("--scores", help="Combined scores", type=str)
    parser.add_argument("--window_size", "-w", help="Window size", type=int, default=50)
    parser.add_argument(
        "--plotfiles_all",
        help="Normalized plot files for all conditions (original.plot.gz).",
        nargs='+',
        type=str,
    )

    # New Algo Params
    parser.add_argument("--disable_new_algorithm", "-disable_newalgo", 
                    help="Disables new categorization and dynamic window for 3,5 prime end generation (default: False)", 
                    action='store_true')
    
    parser.add_argument(
        "--forward_count_condition",
        help="Insertion count tsv file for all conditions (forward.count.tsv).",
        nargs='+',  # Accept one or more arguments as a list
        type=str,
    )
    parser.add_argument(
        "--reverse_count_condition",
        help="Insertion count tsv file for all conditions (reverse.count.tsv).",
        nargs='+',  # Accept one or more arguments as a list
        type=str,
    )
    parser.add_argument(
        "--forward_count_control",
        help="Insertion count tsv file for all controls (forward.count.tsv).",
        nargs='+',  # Accept one or more arguments as a list
        type=str,
    )
    parser.add_argument(
        "--reverse_count_control",
        help="Insertion count tsv file for all controls (reverse.count.tsv).",
        nargs='+',  # Accept one or more arguments as a list
        type=str,
    )
    parser.add_argument(
        "--combined_count_condition",
        help="Insertion count tsv file for all controls (combined.count.tsv).",
        nargs='+',  # Accept one or more arguments as a list
        type=str,
    )
    parser.add_argument("--embl", help="Prepared EMBL file used for analysis", type=str)
    parser.add_argument("--combined_compare", help="Combined compare csvs from comparison folder with logfc, p and q values.", type=str)
    parser.add_argument("--forward_compare", help="Forward compare csv from comparison folder with logfc, p and q values.", type=str)
    parser.add_argument("--reverse_compare", help="Reverse compare csv from comparison folder with logfc, p and q values.", type=str)
    parser.add_argument(
        "--output_dir", "-o", help="Output filename prefix", type=str, default="."
    )
    parser.add_argument(
        "--use_annotations", "-a", help="EMBL file is a real EMBL file with annotations (not a set of windows)", action="store_true", default=False
    )

    
    # Add Gene_Report_Categorization_Params dynamically
    # for param, (arg_name, param_type, default_value) in GENE_REPORT_PARAMS.items():
    #     parser.add_argument(
    #         arg_name,
    #         help=f"{param.replace('_', ' ').capitalize()} (default: {default_value})",
    #         type=param_type
    #     )

    for param, (arg_name, param_type, default_value) in GENE_REPORT_PARAMS.items():
        description = GENE_REPORT_HELP.get(param, param.replace('_', ' ').capitalize())
        parser.add_argument(
            arg_name,
            type=param_type,
            default=default_value,
            help=f"{description} (default: {default_value})"
        )



def gene_report(args):
    """
    Wrapper function to generate a gene essentiality report by calling either the legacy
    or the new version of the `gene_statistics` function based on the `disable_new_algorithm` flag.

    This function integrates normalized insertion count data, comparison statistics,
    gene annotations, and EMBL sequence information to categorize genes based on 
    essentiality, change points, and insertion metrics.

    Args:
        args (argparse.Namespace): Parsed command-line arguments with the following attributes:

            Common:
            - disable_new_algorithm (bool): If True, use the old gene statistics algorithm.
            - embl (str): Path to the prepared EMBL file for gene annotations.
            - output_dir (str): Path to the output directory for saving the gene report.
            - annotations (str, optional): Optional annotation file to override EMBL.
            - use_annotations (bool, optional): If True, use the annotation file instead of the EMBL file.

            Old algorithm only:
            - combined (str): Combined strand plot file.
            - forward (str): Forward strand plot file.
            - reverse (str): Reverse strand plot file.
            - scores (str): Logfc score file, combined.pqvals.plot .
            - window_size (int): Window size used for changepoint detection.

            New algorithm only:
            - plotfiles_all (list of str): List of normalized insertion count plot files for all conditions.
            - forward_count_condition (str): Forward insertion count TSV for the condition.
            - reverse_count_condition (str): Reverse insertion count TSV for the condition.
            - combined_count_condition (str): Combined insertion count TSV for the condition.
            - forward_count_control (str): Forward insertion count TSV for the control.
            - reverse_count_control (str): Reverse insertion count TSV for the control.
            - combined_compare (str): Combined comparison CSV with statistical metrics.
            - forward_compare (str): Forward comparison CSV with statistical metrics.
            - reverse_compare (str): Reverse comparison CSV with statistical metrics.
            - blue_red_logfc_diff_threshold (float): Threshold for logFC difference in blue vs red region.
            - distance_threshold (int): Minimum base-pair distance for changepoint selection.
            - insertion_count_max_threshold (int): Max insertion count threshold.
            - insertion_count_sum_threshold (int): Sum insertion count threshold.
            - insertion_signal_similarity_avg_threshold (float): Average similarity threshold for insertion signal.
            - log_fc_threshold (float): Log fold change threshold.
            - overlap_threshold (float): Overlap threshold for feature regions.
            - q_val (float): Q-value cutoff for significance.
            - inactivation_fraction_threshold (float): Fraction threshold for inactivation determination.

    Returns:
        None. The function writes a gene report to the specified output directory.
    """
    

    if args.disable_new_algorithm:
        gene_statistics_old(
            args.disable_new_algorithm,
            combined_plotfile=args.combined,
            forward_plotfile=args.forward,
            reverse_plotfile=args.reverse,
            combined_scorefile=args.scores,
            window_size=args.window_size,
            embl_file=args.embl,
            output_dir=args.output_dir,
            use_annotations=args.use_annotations,
        )
    else:

        gene_params_keys = [
            "blue_red_logfc_diff_threshold",
            "distance_threshold",
            "insertion_count_max_threshold",
            "insertion_count_sum_threshold",
            "insertion_signal_similarity_avg_threshold",
            "log_fc_threshold",
            "overlap_threshold",
            "q_val",
            "inactivation_fraction_threshold",
        ]
        gene_categorization_params_values = {key: getattr(args, key) for key in gene_params_keys}
        gene_statistics_new(
            args.disable_new_algorithm,
            plotfiles_all=args.plotfiles_all,
            forward_count_condition=args.forward_count_condition,
            reverse_count_condition=args.reverse_count_condition,
            combined_count_condition=args.combined_count_condition,
            forward_count_control=args.forward_count_control,
            reverse_count_control=args.reverse_count_control,
            combined_compare_csv=args.combined_compare,
            forward_compare_csv=args.forward_compare,
            reverse_compare_csv=args.reverse_compare,
            embl_file=args.embl,
            output_dir=args.output_dir,
            gene_categorization_params_values=gene_categorization_params_values,
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
        help="Normalized plot files of both condition and control replicates(original.plot.gz).",
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
        default=198,
    )
    parser.add_argument(
        "--window_interval",
        "-l",
        help="Window interval (default: 25)",
        type=int,
        default=25,
    )
     # New Algo Params
    parser.add_argument("--disable_new_algorithm", "-disable_newalgo", 
                    help="Disables new categorization and dynamic window for 3,5 prime end generation (default: False)", 
                    action='store_true')
    parser.add_argument(
        "--window_size", "-w", help="Window size (default: 100)", type=int, default=100
    )
    # Modification 2
    parser.add_argument("--dynamic_window", "-dw", help="Flag to utilize dynamic window algorithm to generate 3,5 Prime_Features (default: True)", action='store_true')

    # for param, (arg_name, param_type, default_value) in DYNAMIC_WINDOW_PARAMS.items():
    #     if param!="dynamic_window":
    #         parser.add_argument(
    #             arg_name,
    #             help=f"{param.replace('_', ' ').capitalize()} (default: {default_value})",
    #             type=param_type
    #         )
    for param, (arg_name, param_type, default_value) in DYNAMIC_WINDOW_PARAMS.items():
        if param != "dynamic_window":
            description = DYNAMIC_WINDOW_HELP.get(param, param.replace('_', ' ').capitalize())
            parser.add_argument(
                arg_name,
                type=param_type,
                default=default_value,
                help=f"{description} (default: {default_value})"
            )

def prepare_embl(args):

    dynamic_params_keys = [
        "drop_ratio_threshold",
        "gap_threshold",
        "initial_win",
        "initial_win_sum_thres",
        "max_window",
        "min_window",
        "moving_average",
        # "prime_feature_size",
    ]

    dynamic_params_values = {key: getattr(args, key, None) for key in dynamic_params_keys}

    pef = PrepareEMBLFile(
        args.plotfile,
        args.minimum_threshold,
        args.window_size,
        args.window_interval,
        args.prime_feature_size,
        args.emblfile,
        args.dynamic_window,
        args.disable_new_algorithm,
        dynamic_params_values
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

