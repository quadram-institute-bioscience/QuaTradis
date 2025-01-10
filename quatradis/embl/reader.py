import re
from Bio import SeqIO
from collections import defaultdict


class EMBLReader:
    def __init__(self, filename):
        self.filename = filename
        self.features_to_ignore = ['source', 'gene']
        self.genome_length = 0
        self.features = self.read_annotation_features()
        self.genes_to_features = self.gene_names_to_features()

    def read_annotation_features(self):
        self.record = SeqIO.read(self.filename, "embl")
        self.genome_length = len(self.record.seq)

        return [f for f in self.record.features if f.type not in self.features_to_ignore]

    def gene_names_to_features(self):
        genes_to_features = {}
        gene_count = defaultdict(int)
        gene_suffix = defaultdict(int)
        for f in self.features:
            gene_name = self.feature_to_gene_name(f)
            gene_count[gene_name] += 1
            suffix = gene_suffix[gene_name]

            if gene_count[gene_name] > 1:
                gene_suffix[gene_name] += 1
                processed_gene_name = f"{gene_name}_{suffix}"
            else:
                processed_gene_name = gene_name

            genes_to_features[processed_gene_name] = f

        return genes_to_features

    def feature_to_gene_name(self, feature):
        gene_name_val = str(feature.location.start) + "_" + str(feature.location.end)
        if "gene" in feature.qualifiers:
            gene_name_val = feature.qualifiers["gene"][0]
        return gene_name_val

    def calculate_overlap_percentage(self, start1, end1, start2, end2):
        """
        Calculate the percentage of overlap between two ranges.
        """
        overlap_start = max(start1, start2)
        overlap_end = min(end1, end2)
        if overlap_start <= overlap_end:  # Valid overlap
            overlap_length = overlap_end - overlap_start + 1
            # total_length = end1 - start1 + 1
            total_length = end2 - start2 + 1
            return (overlap_length / total_length)
        return 0

    def find_overlaps(self, gene_name):
        """
        Find genes overlapping with the given gene and their overlap percentages.
        Skip genes with __3prime or __5prime in their names.
        """
        if gene_name not in self.genes_to_features:
            raise ValueError(f"Gene {gene_name} not found in the EMBL file.")

        target_feature = self.genes_to_features[gene_name]
        target_start = target_feature.location.start
        target_end = target_feature.location.end

        overlaps = []
        res_given_gene=re.search("^(.+)__([35])prime$", gene_name)

        for other_gene_name, other_feature in self.genes_to_features.items():
            if other_gene_name == gene_name:
                continue  # Skip the same gene

            # Skip genes with __3prime or __5prime in their names
            res = re.search("^(.+)__([35])prime$", other_gene_name)
            if res or res_given_gene.group(1)==other_gene_name:
                continue

            other_start = other_feature.location.start
            other_end = other_feature.location.end
           
            overlap_percentage = self.calculate_overlap_percentage(
                target_start, target_end, other_start, other_end
            )

            if overlap_percentage>0:
                overlaps.append((other_gene_name, overlap_percentage))

        return overlaps
