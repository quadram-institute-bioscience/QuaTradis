"""
Combining insert site plot files
"""
import csv
import os

from Bio import bgzf

from quatradis.tisp.parser import PlotFiles


def stats_header():
    return ["ID", "Sequence Length", "Unique Insertion Sites", "Seq Len/UIS"]


def get_plot_ids(plot_list_file):
    id_order = []
    with open(plot_list_file, 'r') as plf_handle:
        tsv_file = csv.reader(plf_handle, delimiter="\t")
        for row in tsv_file:
            id_order.append(row[0])
    return id_order


def prepare_and_create_tabix_for_combined_plots(tabix_plot, combined_dir):

    tabix_plot_name = combined_dir + os.sep + "tabix.insert_site_plot.gz"
    sorted_tabix_plot_name = combined_dir + os.sep + "tabix_sorted.insert_site_plot.gz"

    with bgzf.BgzfWriter(tabix_plot_name, "wb") as tabix_plot_fh:
        for line in tabix_plot:
            tabix_plot_fh.write(line + "\n")

    os.system("zcat " + tabix_plot_name + " | sort -k1,1 -k2,2n | bgzip > " + sorted_tabix_plot_name +
              " && tabix -b 2 -e 2 " + sorted_tabix_plot_name)
    os.remove(tabix_plot_name)


def combine(plot_file_list, combined_dir="combined"):

    with open(plot_file_list, 'r') as plf_handle:

        stats_file = os.path.splitext(os.path.basename(plot_file_list))[0] + ".stats"
        with open(stats_file, 'w') as stats_file_handle:

            os.makedirs(combined_dir, exist_ok=True)

            stats_writer = csv.writer(stats_file_handle, lineterminator="\n")
            stats_writer.writerow(stats_header())

            plf_reader = csv.reader(plf_handle, delimiter="\t")

            tabix_plot = []
            for line in plf_reader:
                pf = PlotFiles(line[0], line[1:], os.path.dirname(plot_file_list))
                tabix_data, length, nb_uis = pf.combine_plot_files(combined_dir)
                tabix_plot.extend(tabix_data)

                sl_per_uis = str(length / nb_uis) if nb_uis > 0 else "NaN"
                stats_writer.writerow([line[0], length, nb_uis, sl_per_uis])

            if len(tabix_plot) > 0:
                prepare_and_create_tabix_for_combined_plots(tabix_plot, combined_dir)
