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
      extras_require={
        "dev": ["semantic-version", "pytest-cov"],
        "test": ["pytest-cov"]
      },
      scripts=['scripts/tradis',
               'scripts/tradis_comparison.R',
               'scripts/tradis_essentiality.R',
               'pipelines/multi_tradis.nf',
               'pipelines/nextflow.config'])