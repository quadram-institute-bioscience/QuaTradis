import numpy

from Bio.SeqFeature import SeqFeature, FeatureLocation

class Gene:
    def __init__(self, feature=None, blocks=[]):
        self.feature = feature
        self.blocks = blocks
        self.categories = []
        self.upstream = []
        self.expression = ""
        self.direction = ""
        self.five_prime = None
        self.three_prime = None
        self.max_logfc = 0.0
        self.min_pvalue = 1.0
        self.min_qvalue = 1.0
        self.gene_name = self.calc_gene_name() if not self.feature is None else ""

        # Used for merging windows
        self.start = feature.location.start if feature and feature.location else 0
        self.end = feature.location.end if feature and feature.location else 0

    @staticmethod
    def parse_line(line):
        parts = line.split('\t')
        g = Gene()
        g.gene_name = parts[0]
        g.categories = [parts[1]]
        g.feature = SeqFeature(location=FeatureLocation(int(parts[2]), int(parts[3])))
        g.max_logfc = float(parts[4])
        g.expression = parts[5]
        g.direction = parts[6]
        if len(parts) > 7:
            g.upstream = [parts[7]]
            if len(parts) > 8:
                g.min_pvalue = float(parts[8])
                g.min_qvalue = float(parts[9])

        return g

    def calc_gene_name(self):

        if not self.feature:
            return "unknown"

        gene_name_val = str(self.feature.location.start) + "_" + str(self.feature.location.end)
        if "gene" in self.feature.qualifiers:
            gene_name_val = self.feature.qualifiers["gene"][0]
        return gene_name_val


    def category(self):
        if 'knockout' in self.categories:
            return 'knockout'
        elif 'upregulated' in self.categories:
            return 'upregulated'
        elif 'downregulated' in self.categories:
            return 'downregulated'
        else:
            return "/".join(list(set(self.categories)))

    def upstream_gene(self):
        return "/".join(list(set(self.upstream)))

    def max_logfc_from_blocks(self):
        if self.blocks:
            all_logfc = [b.max_logfc for b in self.blocks]
            highest_logfc = numpy.max(numpy.absolute(all_logfc))
            for a in all_logfc:
                if a == highest_logfc:
                    return highest_logfc
                elif a == highest_logfc * -1.0:
                    return highest_logfc * -1.0
        else:
            return 0

    def min_pvalue_from_blocks(self):
        if self.blocks:
            all_pvals = [b.min_pvalue for b in self.blocks]
            return numpy.min(all_pvals)
        else:
            return 0.0

    def min_qvalue_from_blocks(self):
        if self.blocks:
            all_qvals = [b.min_qvalue for b in self.blocks]
            return numpy.min(all_qvals)
        else:
            return 0.0

    def max_logfc_from_category(self):
        max_logfc = self.max_logfc_from_blocks()

        if (max_logfc < 0 and (self.category() == 'upregulated' or self.expression_from_blocks() == 'increased_insertions')) or \
                (max_logfc > 0 and (self.category() == 'downregulated' or self.expression_from_blocks() == 'decreased_insertions')):
            return max_logfc * -1

        return max_logfc

    def expression_from_blocks(self):
        if self.blocks:
            return "/".join(list(set([b.expression for b in self.blocks])))
        else:
            return ""

    def direction_from_blocks(self):
        if self.blocks:
            return "/".join(list(set([b.direction for b in self.blocks])))
        else:
            return ""

    def update_feature(self):
        self.feature = """FT   CDS             {window_start}..{window_end}
    		FT                   /gene="{gene_name}"
    		FT                   /locus_tag="{gene_name}"
    		FT                   /product=product
    		""".format(gene_name=self.gene_name, window_start=str(self.start), window_end=str(self.end))

    def window_string(self):
        return "\t".join(
            [str(self.start) + '_' + str(self.end), str(self.category()), str(self.start), str(self.end),
             str(self.max_logfc_from_category()), str(self.expression_from_blocks()), str(self.direction_from_blocks()),
             str(self.upstream_gene()), str(self.min_pvalue_from_blocks()), str(self.min_qvalue_from_blocks())])

    def simple_string(self):
        return "\t".join([str(self.gene_name), str(self.category()), str(self.feature.location.start),
                                str(self.feature.location.end), str(self.max_logfc),
                                str(self.expression), str(self.direction),
                                str(self.upstream_gene()), str(self.min_pvalue),
                                str(self.min_qvalue)])

    def report_set_string(self, conflict=False):
        return "\t".join([str(self.gene_name), str(self.category()) if not conflict else "conflict", str(self.feature.location.start),
                          str(self.feature.location.end), str(self.max_logfc),
                          str(self.expression), str(self.direction),
                          str(self.upstream_gene()), str(self.min_pvalue),
                          str(self.min_qvalue)])

    def __str__(self):
        try:
            teststring = "\t".join([str(self.gene_name), str(self.category()), str(self.feature.location.start),
                                        str(self.feature.location.end), str(self.max_logfc_from_category()),
                                        str(self.expression_from_blocks()), str(self.direction_from_blocks()),
                                        str(self.upstream_gene()), str(self.min_pvalue_from_blocks()),
                                        str(self.min_qvalue_from_blocks())])
            return teststring
        except:
            raise ValueError("some problem")

    @staticmethod
    def header():
        return "\t".join(['Gene', 'Category', 'Start', 'End', 'MaxLogFC', 'Expression', 'Direction', 'Upstream', 'P-Value', 'Q-Value'])

    @staticmethod
    def report_set_header():
        return "\t".join(['Gene', 'Category', 'Start', 'End', 'Direction', 'Upstream'])
