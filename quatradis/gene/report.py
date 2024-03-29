from quatradis.gene.gene import Gene

class GeneReport:
    '''Read in a Gene report CSV file and return the logFC values against a set list of genes'''

    def __init__(self, filename):
        self.filename = filename
        self.gene_all_data = self.read_all_data_csv()
        self.gene_data = self.read_csv(self.gene_all_data)

    def read_all_data_csv(self):
        with open(self.filename) as genereportfile:
            all_data = [Gene.parse_line(r.strip()) for r in genereportfile if not r.strip().startswith('Gene')]
        return all_data

    def fix_sign_on_logfc(self, gene_all_data):
        for r in gene_all_data:
            if (r.max_logfc < 0.0 and (r.categories[0] == 'upregulated' or r.expression == 'increased_insertions')) or \
                    (r.max_logfc > 0.0 and (
                            r.categories[0] == 'downregulated' or r.expression == 'decreased_insertions')):
                r.max_logfc *= -1.0
        return gene_all_data

    def read_csv(self, gene_all_data):
        return {r.gene_name: r for r in self.fix_sign_on_logfc(gene_all_data)}

    def filtered_genes(self, gene_names):
        row = []
        for gene_name in gene_names:
            if gene_name in self.gene_data:
                row.append(self.gene_data[gene_name])
            else:
                row.append(None)
        return row

    def genes_to_logfc(self, gene_names):
        row = []
        for gene_name in gene_names:
            if gene_name in self.gene_data:
                row.append(str(self.gene_data[gene_name].max_logfc))
            else:
                row.append(str(0.0))
        return row

    def genes_to_qvals(self, gene_names):
        row = []
        for gene_name in gene_names:
            if gene_name in self.gene_all_data:
                row.append(str(self.gene_all_data[gene_name].min_qvalue))
            else:
                row.append(str(1.0))
        return row
