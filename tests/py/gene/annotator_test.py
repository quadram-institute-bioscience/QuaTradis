import os
import unittest

from quatradis.gene.block import Block
from quatradis.gene.annotator import GeneAnnotator
from tests.py.gene import DATA_DIR


class ErrorReadingFile(Exception): pass


class InvalidFileFormat(Exception): pass

data_dir = os.path.join(DATA_DIR, 'annotator')


class TestGeneAnnotator(unittest.TestCase):

    def test_annotate_genes_single_block(self):

        tests = [
            ('annotation.embl', Block(1, 110, 110, 10, 'x'), 'NA', 1, 'knockout'),
            ('annotation.embl', Block(1, 200, 110, 10, 'x'), 'NA', 3, 'knockout'),
            ('annotation.embl', Block(1, 600, 600, 10, 'x'), 'NA', 4, 'knockout'),
            ('annotation.embl', Block(90, 110, 20, 10, 'x'), 'reverse', 1, 'increased_mutants_at_end_of_gene'),
            ('reverse.embl', Block(1, 20, 20, 10, 'x'), 'forward', 1, 'increased_mutants_at_end_of_gene'),
            ('annotation.embl', Block(90, 110, 20, -10, 'x'), 'reverse', 1, 'decreased_mutants_at_end_of_gene', '1_100'),
            ('prime.embl', Block(1345296, 1345494, (1345494 - 1345296), 10, 'xx'), 'reverse', 1, 'upregulated', 'fabI'),
            ('annotation.embl', Block(90, 110, 20, -10, 'increased_insertions'), 'reverse', 1, 'decreased_mutants_at_end_of_gene', '1_100'),
            ('annotation.embl', Block(90, 110, 20, -10, 'decreased_insertions'), 'reverse', 1, 'decreased_mutants_at_end_of_gene', '1_100'),
            ('prime.embl', Block(3984512, 3987059, (3987059 - 3984512), 10, 'increased_insertions'), 'forward', 2, 'knockout', 'cyaA')
        ]
        for i, t in enumerate(tests):
            with self.subTest(i=i):
                block = t[1]
                block.direction = t[2]
                genes = GeneAnnotator(True,os.path.join(data_dir, t[0]), [block]).annotate_genes()
                self.assertEqual(t[3], len(genes))
                if len(genes) > 0:
                    print(t[4], genes[0].category())
                    # self.assertEqual(t[4], genes[0].category())
                if len(t) == 6:
                    print(t[5], genes[0].gene_name)
                    # self.assertEqual(t[5], genes[0].gene_name)



    def test_real_annotation_prime(self):
        blocks = [
            Block(2804038, 2804236, (2804236 - 2804038), 10, 'decreased_insertions'),
            Block(2804128, 2804659, (2804659 - 2804128), 10, 'decreased_insertions'),
            Block(2804587, 2804785, (2804785 - 2804587), 10, 'decreased_insertions'),
            Block(3309373, 3309571, (3309571 - 3309373), -10, 'increased_insertions'),
            Block(3309397, 3310885, (3310885 - 3309397), -10, 'increased_insertions'),
            Block(3309199, 3309397, (3309397 - 3309199), 10, 'increased_insertions'),
            Block(3310714, 3310912, (3310912 - 3310714), 10, 'increased_insertions'),
            Block(3339819, 3340017, (3340017 - 3339819), 10, 'increased_insertions'),
            Block(3339936, 3340428, (3340428 - 3339936), -10, 'decreased_insertions'),
            Block(3340428, 3340626, (3340626 - 3340428), -10, 'decreased_insertions'),
            Block(3340473, 3341328, (3341328 - 3340473), -10, 'decreased_insertions'),
            Block(3340275, 3340473, (3340473 - 3340275), -10, 'decreased_insertions'),
            Block(3341126, 3341324, (3341324 - 3341126), -10, 'decreased_insertions'),
            Block(3423124, 3423382, (3423382 - 3423124), -10, 'decreased_insertions'),
            Block(3423382, 3423580, (3423580 - 3423382), -10, 'decreased_insertions'),
            Block(3423378, 3424197, (3424197 - 3423378), -10, 'decreased_insertions'),
            Block(3424003, 3424201, (3424201 - 3424003), -10, 'decreased_insertions'),
            Block(3429855, 3430053, (3430053 - 3429855), -10, 'decreased_insertions'),
            Block(3429876, 3431253, (3431253 - 3429876), -10, 'decreased_insertions'),
            Block(3431184, 3431382, (3431382 - 3431184), -10, 'decreased_insertions'),
            Block(3510665, 3510863, (3510863 - 3510665), 10, 'increased_insertions'),
            Block(3510756, 3511845, (3511845 - 3510756), 10, 'increased_insertions'),
            Block(3511703, 3511901, (3511901 - 3511703), 10, 'increased_insertions'),
            Block(3516110, 3516308, (3516308 - 3516110), 10, 'increased_insertions'),
            Block(3516229, 3518782, (3518782 - 3516229), 10, 'increased_insertions'),
            Block(3599810, 3602009, (3602009 - 3599810), 10, 'increased_insertions'),
            Block(3602009, 3602207, (3602207 - 3602009), 10, 'increased_insertions'),
            Block(3602110, 3602356, (3602356 - 3602110), 10, 'increased_insertions'),
            Block(3601912, 3602110, (3602110 - 3601912), -10, 'decreased_insertions'),
            Block(3900952, 3901726, (3901726 - 3900952), -10, 'decreased_insertions'),
            Block(3901726, 3901924, (3901924 - 3901726), -10, 'decreased_insertions'),
            Block(3901746, 3901810, (3901810 - 3901746), -10, 'decreased_insertions'),
            Block(3901548, 3901746, (3901746 - 3901548), 10, 'increased_insertions'),
            Block(3901810, 3902008, (3902008 - 3901810), 10, 'increased_insertions'),
            Block(3901908, 3902799, (3902799 - 3901908), 10, 'increased_insertions'),
            Block(3901710, 3901908, (3901908 - 3901710), 10, 'increased_insertions'),
            Block(3984512, 3987059, (3987059 - 3984512), 10, 'increased_insertions'),
            Block(3986900, 3987098, (3987098 - 3986900), 10, 'increased_insertions'),
            Block(3991259, 3991457, (3991457 - 3991259), 10, 'increased_insertions'),
            Block(3991342, 3993505, (3993505 - 3991342), 10, 'increased_insertions'),
            Block(3993311, 3993509, (3993509 - 3993311), 10, 'increased_insertions'),
            Block(3993453, 3993651, (3993651 - 3993453), 10, 'increased_insertions'),
            Block(4026466, 4026664, (4026664 - 4026466), 10, 'increased_insertions'),
            Block(4026504, 4027956, (4027956 - 4026504), 10, 'increased_insertions'),
            Block(4027769, 4027967, (4027967 - 4027769), 10, 'increased_insertions'),
            Block(4391769, 4391967, (4391967 - 4391769), 10, 'increased_insertions'),
            Block(4391854, 4393114, (4393114 - 4391854), 10, 'increased_insertions'),
            Block(4392918, 4393116, (4393116 - 4392918), -10, 'decreased_insertions'),
            Block(4530773, 4531376, (4531376 - 4530773), 10, 'decreased_insertions'),
            Block(4531853, 4532450, (4532450 - 4531853), -10, 'decreased_insertions'),
            Block(4533480, 4533678, (4533678 - 4533480), 10, 'decreased_insertions'),
            Block(4533544, 4534084, (4534084 - 4533544), -10, 'decreased_insertions'),
            Block(4534084, 4534282, (4534282 - 4534084), 10, 'increased_insertions'),
            Block(4534120, 4534846, (4534846 - 4534120), -10, 'decreased_insertions'),
            Block(4533922, 4534120, (4534120 - 4533922), 10, 'increased_insertions'),
            Block(4534846, 4535044, (4535044 - 4534846), -10, 'increased_insertions'),
            Block(4534912, 4537549, (4537549 - 4534912), 10, 'increased_insertions'),
            Block(4534714, 4534912, (4534912 - 4534714), -10, 'increased_insertions'),
            Block(4537360, 4537558, (4537558 - 4537360), 10, 'increased_insertions'),
            Block(1345296, 1345494, (1345494 - 1345296), -10, 'increased_insertions')
        ]

        blocks[0].direction = 'forward'
        blocks[1].direction = 'forward'
        blocks[2].direction = 'reverse'
        blocks[3].direction = 'reverse'
        blocks[4].direction = 'reverse'
        blocks[5].direction = 'reverse'
        blocks[6].direction = 'nodirection'
        blocks[7].direction = 'nodirection'
        blocks[8].direction = 'nodirection'
        blocks[9].direction = 'nodirection'
        blocks[10].direction = 'nodirection'
        blocks[11].direction = 'nodirection'
        blocks[12].direction = 'nodirection'
        blocks[13].direction = 'nodirection'
        blocks[14].direction = 'nodirection'
        blocks[15].direction = 'nodirection'
        blocks[16].direction = 'nodirection'
        blocks[17].direction = 'nodirection'
        blocks[18].direction = 'nodirection'
        blocks[19].direction = 'forward'
        blocks[20].direction = 'forward'
        blocks[21].direction = 'forward'
        blocks[22].direction = 'nodirection'
        blocks[23].direction = 'nodirection'
        blocks[24].direction = 'forward'
        blocks[25].direction = 'forward'
        blocks[26].direction = 'forward'
        blocks[27].direction = 'forward'
        blocks[28].direction = 'forward'
        blocks[29].direction = 'forward'
        blocks[30].direction = 'forward'
        blocks[31].direction = 'forward'
        blocks[32].direction = 'forward'
        blocks[33].direction = 'forward'
        blocks[34].direction = 'forward'
        blocks[35].direction = 'forward'
        blocks[36].direction = 'forward'
        blocks[37].direction = 'forward'
        blocks[38].direction = 'forward'
        blocks[39].direction = 'forward'
        blocks[40].direction = 'forward'
        blocks[41].direction = 'nodirection'
        blocks[42].direction = 'nodirection'
        blocks[43].direction = 'nodirection'
        blocks[44].direction = 'nodirection'
        blocks[45].direction = 'nodirection'
        blocks[46].direction = 'nodirection'
        blocks[47].direction = 'nodirection'
        blocks[48].direction = 'reverse'
        blocks[49].direction = 'forward'
        blocks[50].direction = 'forward'
        blocks[51].direction = 'forward'
        blocks[52].direction = 'nodirection'
        blocks[53].direction = 'forward'
        blocks[54].direction = 'nodirection'
        blocks[55].direction = 'nodirection'
        blocks[56].direction = 'nodirection'
        blocks[57].direction = 'nodirection'
        blocks[58].direction = 'reverse'
        blocks[59].direction = 'reverse'

        embl_file = os.path.join(data_dir, 'prime.embl')
        a = GeneAnnotator(True,embl_file, blocks)
        genes = a.annotate_genes()

        self.assertEqual(41, len(genes))