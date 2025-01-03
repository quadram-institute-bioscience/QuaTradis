import os
import unittest

from quatradis.artemis.experiment import ExperimentCollection, ExperimentMetaData
from tests.py.artemis import DATA_DIR

data_dir = os.path.join(DATA_DIR, 'experimentcollection')


class TestExperimentCollection(unittest.TestCase):

    def test_vary_mic(self):
        e = [ExperimentMetaData('drug', 'X', 'X', 'X', 1, 0, os.path.join(data_dir, 'condition_1_0_rep1'),
                                os.path.join(data_dir, 'condition_1_0_rep2')),
             ExperimentMetaData('drug', 'X', 'X', 'X', 0.5, 0, os.path.join(data_dir, 'condition_0.5_0_rep1'),
                                os.path.join(data_dir, 'condition_0.5_0_rep2')),
             ExperimentMetaData('drug', 'X', 'X', 'X', 2, 0, os.path.join(data_dir, 'condition_2_0_rep1'),
                                os.path.join(data_dir, 'condition_2_0_rep2')),
             ExperimentMetaData('drug', 'X', 'X', 'X', 0.25, 0, os.path.join(data_dir, 'condition_0.25_0_rep1'),
                                os.path.join(data_dir, 'condition_0.25_0_rep2'))]
        ec = ExperimentCollection(e, [os.path.join(data_dir, 'control')], os.path.join(data_dir, 'reference'))
        self.assertEqual({'drug': 'drug', 'induction': 0}, ec.properties_in_common())
        self.assertEqual(['mic'], ec.properties_not_in_common())
        self.assertEqual('drug_0', ec.project_name)
