import os

from setuptools import setup

requirements = [line.rstrip() for line in open("requirements.txt", "rt")]


def read(fname):
    return open(os.path.join(os.path.dirname(__file__), fname)).read()


setup(name='quatradis',
      version=read('VERSION'),
      description='Analyses transposon insertion sites',
      long_description=read('README.md'),
      url='https://quadram.ac.uk',
      author='Sarah Bastkowski',
      author_email='sarah.bastkowski@quadram.ac.uk',
      license=read('LICENSE'),
      packages=['quatradis'],
      zip_safe=False,
      install_requires=requirements,
      scripts=[   'scripts/tradis',
                  'scripts/tradis_nf',
                  'scripts/tradis_plot',
                  'scripts/combine_tradis_plots',
                  'scripts/add_tradis_tags',
                  'scripts/check_tradis_tags',
                  'scripts/remove_tradis_tags',
                  'scripts/tradis_gene_insert_sites',
                  'scripts/index_reference',
                  'scripts/tradis_comparison.R',
                  'scripts/tradis_essentiality.R'])
