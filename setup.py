import os

from setuptools import setup, Extension
from Cython.Build import cythonize

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
      ext_modules=cythonize([Extension("quatradis.tisp.write", ["quatradis/tisp/write.pyx"])]),
      zip_safe=False,
      install_requires=requirements,
      extras_require={
        "dev": ["semantic-version", "pytest-cov"],
        "test": ["pytest-cov"]
      },
      scripts=['tradis',
               'quatradis/comparison/tradis_comparison.R',
               'quatradis/essentiality/tradis_essentiality.R',
               'quatradis/pipelines/create_plots.smk',
               'quatradis/pipelines/compare.smk'])