import unittest
import os

from quatradis.gene.report import GeneReport


class ErrorReadingFile(Exception): pass


class InvalidFileFormat(Exception): pass


data_dir = os.path.join('data', 'gene', 'report')


class TestGeneReport(unittest.TestCase):

    def test_valid_file(self):
        g = GeneReport(os.path.join(data_dir, 'ctrl.csv'))
        self.assertTrue(g)
        gene_names = ['perR', 'ykfC', 'ykgH', 'yaiL', 'yaiP', 'unknown']
        self.assertEqual(['10.0', '11.0', '10.0', '10.0', '11.0', '10.0'], g.genes_to_logfc(gene_names))

        # reorder the gene names
        gene_names = ['yaiL', 'ykfC', 'ykgH', 'perR', 'unknown', 'yaiP']
        self.assertEqual(['10.0', '11.0', '10.0', '10.0', '10.0', '11.0'], g.genes_to_logfc(gene_names))
