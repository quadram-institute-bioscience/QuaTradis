import os

from quatradis.gene.report import GeneReport
from quatradis.gene.gene import Gene


def gene_reports_run(args):
    if not os.path.exists(args.prefix):
        os.makedirs(args.prefix)
    p = GeneReportSets(args.genereports, args.prefix)
    p.write_union_file()


class GeneReportSets:
    '''Take in 2 or more gene report spreadsheets and output the union and intersection, ...'''

    def __init__(self, filenames, prefix):
        self.filenames = filenames
        self.prefix = prefix

        if not os.path.exists(self.prefix):
            os.makedirs(self.prefix)

        self.gene_reports = self.parse_gene_reports()

    def parse_gene_reports(self):
        return {filename: GeneReport(filename) for filename in self.filenames}

    # Find the genes that are present in each file
    def intersection(self):

        # Not super efficient but neater and should be fast enough
        gene_sets = []
        gene_category_sets = []
        combined, combined_with_categories = self.union()
        for i, f in enumerate(self.filenames):
            genes = set()
            gene_categories = set()
            for gene in self.gene_reports[f].gene_all_data:
                genes.add(gene.gene_name)
                gene_categories.add(gene.gene_name + '~' + gene.category())

            gene_sets.append(genes)
            gene_category_sets.append(gene_categories)
            print(f, ": Genes=", len(genes), "; Genes+Categories=", len(gene_categories))

        gene_intersection = set.intersection(*gene_sets)
        gene_category_intersection = set.intersection(*gene_category_sets)

        print("Gene Intersection =", len(gene_intersection))
        print("Gene+Category Intersection =", len(gene_category_intersection))

        filtered_combined = {}
        filtered_categories = {}

        for gene in gene_intersection:
            filtered_combined[gene] = combined[gene]

        for gene_category in gene_category_intersection:
            if gene_category not in combined_with_categories:
                print("Couldn't find", gene_category, "in union.")
            else:
                filtered_categories[gene_category] = combined_with_categories[gene_category]

        return filtered_combined, filtered_categories

    @staticmethod
    def row_to_gene_name(row):
        gene = row[0]
        # use the start and end coords for unnamed features
        if gene == 'unknown' or gene == 'NA':
            gene = str(row[2]) + "_" + str(row[3])
        return gene

    # When merging, use the first row for a gene. This can have unintended consequences (like an increase in insertions in one exp, and a decrease in insertions in another)
    def union(self):
        combined = {}
        combined_with_categories = {}
        for f in self.filenames:
            for gene in self.gene_reports[f].gene_all_data:
                if gene.gene_name not in combined:
                    combined[gene.gene_name] = gene
                combined_with_categories[gene.gene_name + '~' + gene.category()] = gene

        return combined, combined_with_categories

    def write_union_file(self):
        union_filename = os.path.join(self.prefix, "union_gene_report.csv")

        union, union_categories = self.union()
        is_conflict = []

        for i in union_categories:
            is_conflict.append(i.split("~")[0])

        with open(union_filename, 'w') as bf:
            bf.write(Gene.header() + "\n")
            for gene in sorted(union.values(), key=lambda x: x.feature.location.start):
                num_conflict = is_conflict.count(str(gene.gene_name))
                conflict = num_conflict > 1
                bf.write(str(gene.report_set_string(conflict=conflict)) + "\n")

        return self

    def write_intersection_file(self):
        intersection_filename = os.path.join(self.prefix, "intersection_gene_report.csv")
        intersection, intersection_categories = self.intersection()
        if len(intersection) > 0:
            with open(intersection_filename, 'w') as bf:
                bf.write(str(Gene.header()) + "\n")
                for gene in sorted(intersection.values(), key=lambda x: x.feature.location.start):
                    conflict = str(gene.gene_name + '~' + gene.category()) not in intersection_categories
                    bf.write(str(gene.report_set_string(conflict=conflict)) + "\n")
        else:
            print("No intersecting genes")

        return self
