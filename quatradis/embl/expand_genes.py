""" Given an annotation file, take each gene, and create a new feature at the start and end to capture promotors"""

from quatradis.embl.reader import EMBLReader
from quatradis.embl.sequence import EMBLSequence
import numpy as np
import pandas as pd


class FeatureProperties:
    def __init__(
        self, start=0, end=0, direction="", gene_name="", locus_tag="", product=""
    ):
        self.start = start
        self.end = end
        self.direction = direction
        self.gene_name = gene_name
        self.locus_tag = locus_tag
        self.product = product


class   EMBLExpandGenes:
    def __init__(self, embl_file,prime_feature_size,dynamic_window,kwargs):
        # Modification 6
        self.embl_file = embl_file
        self.dynamic_window=dynamic_window
        # Initialize dynamic parameters here
        self.drop_ratio_threshold = kwargs.get("drop_ratio_threshold", None)
        self.gap_threshold = kwargs.get("gap_threshold", None)
        self.initial_win = kwargs.get("initial_win", None)
        self.initial_win_sum_thres = kwargs.get("initial_win_sum_thres", None)
        self.max_window = kwargs.get("max_window", None)
        self.min_window = kwargs.get("min_window", None)
        self.moving_average = kwargs.get("moving_average", None)
        # if self.dynamic_window:
        #     self.feature_size = kwargs.get("prime_feature_size",None)
        # else:
        self.feature_size=prime_feature_size

        self.er = EMBLReader(self.embl_file)
        self.features = self.er.features
        self.genes_to_features= self.er.genes_to_features
        self.genome_length = self.er.genome_length
        

    def create_3_5_prime_features(self):
        new_features = []
        for gene, feature in self.genes_to_features.items():
            gene_name= gene
            locus_tag = self.feature_to_locus_tag(feature)
            product = self.feature_to_product(feature)

            # The gene itself
            new_features.append(
                FeatureProperties(
                    feature.location.start,
                    feature.location.end,
                    feature.location.strand,
                    gene_name,
                    locus_tag,
                    product,
                )
            )

            # forward direction
            if feature.location.strand == 1:
                new_features.append(
                    self.construct_start_feature(
                        feature, gene_name, "__5prime", locus_tag, product
                    )
                )
                new_features.append(
                    self.construct_end_feature(
                        feature, gene_name, "__3prime", locus_tag, product
                    )
                )
            else:
                new_features.append(
                    self.construct_end_feature(
                        feature, gene_name, "__5prime", locus_tag, product
                    )
                )
                new_features.append(
                    self.construct_start_feature(
                        feature, gene_name, "__3prime", locus_tag, product
                    )
                )

        return new_features
    

    def get_windows(self,prime_type, data, start, end, max_window, min_window):
        #red inserts is zero, blue is one
        low = max(1, start - max_window)
        high = min(data.shape[0], end + max_window) 

        #Indices are preserved when slicing, even after reversing
        if prime_type=="start": 
            shifts = self.get_shifts(data.loc[low:start-1].iloc[:, 0].iloc[::-1])
            if (not shifts.empty) and (shifts[0] < (start - self.initial_win)):
                if(data.loc[(start - self.initial_win):(start - self.initial_win + 99)].iloc[:, 0].sum() < self.initial_win_sum_thres):
                    return max(1, start - min_window)
            
            gap_in_shifts = shifts.diff(periods=1)
            for index in gap_in_shifts.index:
                if gap_in_shifts[index] > self.gap_threshold:
                    return min(min(shifts[0:index]),(max(1, start - min_window)))
            
            if shifts.empty:
                return max(1, start - min_window)
            else:
                return min(min(shifts),(max(1, start - min_window)))

        elif prime_type=="end":
            shifts = self.get_shifts(data.loc[end+1:high].iloc[:, 1])
            if (not shifts.empty) and (shifts[0] > (end + self.initial_win)):
                if(data.loc[(end + self.initial_win - 99):(end + self.initial_win)].iloc[:, 1].sum() < self.initial_win_sum_thres):
                    return min(data.shape[0], end + min_window)
            gap_in_shifts = shifts.diff(periods=1)
            for index in gap_in_shifts.index:
                if gap_in_shifts[index] > self.gap_threshold:
                    return max(max(shifts[0:index]),(min(data.shape[0], end + min_window)))
                
            if shifts.empty:
                return min(data.shape[0], end + min_window)
            else:
                return max(max(shifts),(min(data.shape[0], end + min_window)))
            

    def get_shifts(self, data, window_size=5):
        moving_average = data.rolling(window=window_size, min_periods=1).mean().shift(1)
        drop_ratio = data/moving_average
        potential_shifts = drop_ratio[(drop_ratio < self.drop_ratio_threshold) & (moving_average > self.moving_average)].index.tolist()
        return pd.Series(potential_shifts)
        

    # Modification 9
    def construct_end_feature(self, feature, gene_name, suffix, locus_tag, product):
        # density_check= self.check_prime_feature_density("end",feature.location.end)
        # if density_check and self.dynamic_window:
        if self.dynamic_window:
            end_index_list=[]
            for key, parser in self.plot_file_objs.items():
                if key.startswith("Condition"):
                    index= self.get_windows("end",parser.insert_site_array, feature.location.start, feature.location.end, self.max_window,self.min_window)
                    end_index_list.append(index)
            start= feature.location.end
            end = min(end_index_list)
        else:
            start = feature.location.end
            end = feature.location.end + self.feature_size

        if end > self.genome_length:
            end = self.genome_length

        if start >= end or end - start < 10:
            return None

        return FeatureProperties(
            start,
            end,
            feature.location.strand,
            gene_name + suffix,
            locus_tag + suffix,
            product,
        )

    # Modification 10
    def construct_start_feature(self, feature, gene_name, suffix, locus_tag, product,dynamic_window=False):

        if self.dynamic_window:
            start_index_list=[]
            for key, parser in self.plot_file_objs.items():
                if key.startswith("Condition"):
                    index= self.get_windows("start",parser.insert_site_array, feature.location.start, feature.location.end, self.max_window,self.min_window)
                    start_index_list.append(index)
            start= max(start_index_list)
            end =feature.location.start
            
        else:
            start = feature.location.start - self.feature_size
            end = feature.location.start
        if start < 1:
            start = 1
        if start >= end or end - start < 10:
            return None
        return FeatureProperties(
            start,
            end,
            feature.location.strand,
            gene_name + suffix,
            locus_tag + suffix,
            product,
        )
    # Modification 11
    def construct_file(self, filename,plot_file_objs):
        with open(filename, "w") as emblfile:
            emblfile.write(self.header())
            self.plot_file_objs= plot_file_objs
            for f in self.create_3_5_prime_features():
                if f == None:
                    continue
                if f.direction == 1:
                    emblfile.write(self.construct_feature_forward(f))
                else:
                    emblfile.write(self.construct_feature_reverse(f))

            emblfile.write(EMBLSequence(str(self.er.record.seq)).format())
        return self

    def feature_to_gene_name(self, feature):
        gene_name_val = str(feature.location.start) + "_" + str(feature.location.end)
        if "gene" in feature.qualifiers:
            gene_name_val = feature.qualifiers["gene"][0]
        return gene_name_val

    def feature_to_product(self, feature):
        product_val = str(feature.location.start) + "_" + str(feature.location.end)
        if "product" in feature.qualifiers:
            product_val1 = feature.qualifiers["product"][0]
            product_val = product_val1.replace(",", " and ")
        return product_val

    def feature_to_locus_tag(self, feature):
        locus_tag_val = str(feature.location.start) + "_" + str(feature.location.end)
        if "locus_tag" in feature.qualifiers:
            locus_tag_val = feature.qualifiers["locus_tag"][0]
        return locus_tag_val

    def header(self):
        return """ID   ABC; SV 1; circular; genomic DNA; STD; PRO; {length} BP.
XX
FH   Key             Location/Qualifiers
FH
FT   source          1..{length}
FT                   /organism="Bacteria"
""".format(
            length=str(self.genome_length)
        )

    def construct_feature_forward(self, feature):
        return """FT   CDS             {window_start}..{window_end}
FT                   /gene="{gene_name}"
FT                   /locus_tag="{locus_tag}"
FT                   /product="{product}"
""".format(
            gene_name=feature.gene_name,
            window_start=str(feature.start + 1),
            window_end=str(feature.end),
            locus_tag=feature.locus_tag,
            product=feature.product,
        )

    def construct_feature_reverse(self, feature):
        return """FT   CDS             complement({window_start}..{window_end})
FT                   /gene="{gene_name}"
FT                   /locus_tag="{locus_tag}"
FT                   /product="{product}"
""".format(
            gene_name=feature.gene_name,
            window_start=str(feature.start + 1),
            window_end=str(feature.end),
            locus_tag=feature.locus_tag,
            product=feature.product,
        )
