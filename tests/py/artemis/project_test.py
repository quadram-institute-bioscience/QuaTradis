import os
import unittest

from quatradis.artemis.project import ArtemisProject
from tests.py.artemis import DATA_DIR

data_dir = os.path.join(DATA_DIR, 'experimentcollection')


class TestArtemisProject(unittest.TestCase):

    def test_vary_mic(self):
        outputfile = "project_file"
        ap = ArtemisProject(outputfile, False, os.path.join(data_dir, 'spreadsheet'),
                            os.path.join(data_dir, 'reference'), [os.path.join(data_dir, 'control')])

        self.assertTrue(ap.create_project_file())
        self.assertTrue(os.path.exists(outputfile))

        os.remove(outputfile)
