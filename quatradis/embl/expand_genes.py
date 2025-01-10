""" Given an annotation file, take each gene, and create a new feature at the start and end to capture promotors"""

from quatradis.embl.reader import EMBLReader
from quatradis.embl.sequence import EMBLSequence
import numpy as np


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
    def __init__(self, embl_file, feature_size,dynamic_window=False):
        # Modification 6
        self.embl_file = embl_file
        self.feature_size = feature_size
        self.er = EMBLReader(self.embl_file)
        self.features = self.er.features
        self.genes_to_features= self.er.genes_to_features
        self.genome_length = self.er.genome_length
        self.dynamic_window=dynamic_window
        self.max_window = 2000
        self.prime_density_threshold= 2
        self.shift_continuity_requirement=150


    def create_3_5_prime_features(self):
        new_features = []
        for gene, feature in self.genes_to_features.items():
            gene_name= gene
            # gene_name = self.feature_to_gene_name(feature)
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
    # Modification 7
    def get_windows(self,prime_type,geneName,data, start, end, max_window):
        #red inserts is zero, blue is one
        low = max(1, start - max_window)
        high = min(data.size, end + max_window)
         
        if prime_type=="start": 
            # blue = self.get_shifts(geneName,data[end:high].iloc[:, 1],)
            shifts = self.get_shifts(data[low:start].iloc[:, 0].iloc[::-1])
            shifts.append(max(1, start - 200))
            if geneName in ["acrA"]:
              print(prime_type)
              print("Start: ",geneName)
              print("shifts",shifts)
            return min(shifts)
        elif prime_type=="end":
            shifts= self.get_shifts(data[end:high].iloc[:, 1])
            shifts.append(min(len(data), end + 200))
            if geneName in ["acrA"]:
              print(prime_type)
              print("End: ",geneName)
              print("shifts",shifts)
            return max(shifts)

            

    # Modification 8
    def get_shifts(self,data, window_size=20):
        moving_average = data.rolling(window=window_size).mean()
        deviation = abs(data - moving_average)
        threshold = 2 * data.std()  
        potential_shifts = deviation[deviation > threshold].index.tolist()
        if len(potential_shifts)>1:
            for indx in range(len(potential_shifts) - 1):
                    diff = abs(potential_shifts[indx+1] - potential_shifts[indx])
                    if diff > self.shift_continuity_requirement:
                        return potential_shifts[0:indx]
        return potential_shifts
    
    def check_prime_feature_density(self,prime_type,gene_index):
        condition_sum = 0
        control_sum = 0
        # print("gene_index",gene_index)
        genome_length= self.plot_file_objs["Condition1"].genome_length
        
        for key, parser in self.plot_file_objs.items():
            if prime_type == 'start':
                if gene_index - self.feature_size <0:
                    start_index=0
                else:
                    start_index=gene_index - self.feature_size
                array_slice = parser.insert_site_array.values[start_index:gene_index, 0]
            elif prime_type == 'end':
                if gene_index + self.feature_size>genome_length:
                    end_index=genome_length
                else:
                    end_index=gene_index + self.feature_size
                array_slice = parser.insert_site_array.values[gene_index:end_index, 1]
            else:
                raise ValueError("prime_type must be 'start' or 'end'")

            if key.startswith("Condition"):
                condition_sum += np.sum(array_slice)
            elif key.startswith("Control"):
                control_sum += np.sum(array_slice)

        if control_sum == 0 and condition_sum!=0:
            return True
        elif control_sum == 0 and condition_sum==0:
            return False

        ratio = condition_sum / control_sum
        return ratio > self.prime_density_threshold

    # Modification 9
    def construct_end_feature(self, feature, gene_name, suffix, locus_tag, product):
        density_check= self.check_prime_feature_density("end",feature.location.end)
        if gene_name=="yciW":
            print("construct_end_feature")
            print("gene_name",gene_name)
            print("density_check",density_check)
        if density_check and self.dynamic_window:
            end_index_list=[]
            for key, parser in self.plot_file_objs.items():
                if key.startswith("Condition"):
                    index= self.get_windows("end",gene_name,parser.insert_site_array, feature.location.start, feature.location.end, self.max_window)
                    end_index_list.append(index)
            if gene_name=="acrA":
                print("end_index_list",end_index_list)
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
        density_check= self.check_prime_feature_density("start",feature.location.start)
        if gene_name=="yciW":
            print("construct_start_feature")
            print("gene_name",gene_name)
            print("density_check",density_check)
        if density_check and self.dynamic_window:
            start_index_list=[]
            for key, parser in self.plot_file_objs.items():
                if key.startswith("Condition"):
                    index= self.get_windows("start",gene_name,parser.insert_site_array, feature.location.start, feature.location.end, self.max_window)
                    start_index_list.append(index)
            if gene_name=="acrA":
                print("start_index_list",start_index_list)
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
