DYNAMIC_WINDOW_PARAMS= {
    "dynamic_window": ("--dynamic_window", bool,True),
    "drop_ratio_threshold": ("--drop_ratio_threshold", float, 0.5),
    "gap_threshold": ("--gap_threshold", int, 200),
    "initial_win": ("--initial_win", int, 400),
    "initial_win_sum_thres": ("--initial_win_sum_thres", int, 200),
    "max_window": ("--max_window", int, 2000),
    "min_window": ("--min_window", int, 198),
    "moving_average": ("--moving_average", int, 50),
    # "prime_feature_size": ("--prime_feature_size", int, 198),
    }

GENE_REPORT_PARAMS = {
        "blue_red_logfc_diff_threshold": ("--blue_red_logfc_diff_threshold", float, 0.4),
        "distance_threshold": ("--distance_threshold", float, 0.6),
        "insertion_count_max_threshold": ("--insertion_count_max_threshold", int, 30),
        "insertion_count_sum_threshold": ("--insertion_count_sum_threshold", int, 75),
        "insertion_signal_similarity_avg_threshold": ("--insertion_signal_similarity_avg_threshold", float, 0.5),
        "log_fc_threshold": ("--log_fc_threshold", float, 2.0),
        "overlap_threshold": ("--overlap_threshold", float, 0.6),
        "q_val": ("--q_val", float, 0.05),
        "inactivation_fraction_threshold": ("--inactivation_fraction_threshold", float, 0.2),
    }

GENE_REPORT_HELP = {
    "blue_red_logfc_diff_threshold": "LogFC difference threshold between forward and reverse strands, used to establish precedence between knockout and up/down regulation.",
    "distance_threshold": "Shape-based similarity measure to capture consistency  across conditions (lower value depicts similar replicates).",
    "insertion_count_max_threshold": "Minimum required average of maximum insertions at a base pair across replicates required for gene categorization.",
    "insertion_count_sum_threshold": "Minimum average total insertions across condition replicates required for gene categorization.",
    "insertion_signal_similarity_avg_threshold": "Maximum allowed variation in average of non-zero insertion values across condition replicates.",
    "log_fc_threshold": "Minimum Log fold change (LogFC) required for protection, knockout, and up/down regulation categorization.",
    "overlap_threshold": "Minimum overlap between an adjacent gene and the prime end required for explaining categorization knockout due to up/down regulation and vice versa.",
    "q_val": "Q-value (modified p-value) significance threshold.",
    "inactivation_fraction_threshold": "Minimum percentage length of gene inactivation required for knockout categorization.",
}
 
DYNAMIC_WINDOW_HELP = {
    "max_window": "Maximum window size that can be created",
    "min_window": "Minimum window size that can be created",
    "initial_win": "Initial window size within which the insertion density is checked before dynamic window propagation",
    "initial_win_sum_thres": "Threshold for the sum of insertions within the initial window",
    "gap_threshold": "Maximum allowed discontinuity in the insertion sequence in the window",
    "drop_ratio_threshold": "Percentage drop in insertions required to be considered a change point",
    "moving_average": "Moving average of indices with a window size of 5"
}