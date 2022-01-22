import os

from setuptools import setup

requirements = [line.rstrip() for line in open("requirements.txt", "rt")]


def read(fname):
    return open(os.path.join(os.path.dirname(__file__), fname)).read()


setup(name='quatradis',
      version=read('VERSION'),
      description='Analyses transposon insertion sites',
      long_description=read('README.md'),
      url='http://quadram.ac.uk',
      author='Sarah Bastkowski',
      author_email='sarah.bastkowski@quadram.ac.uk',
      license=read('LICENSE'),
      packages=['quatradis'],
      zip_safe=False,
      install_requires=requirements,
      scripts=[   'scripts/bacteria_tradis',
                  'scripts/tradis_plot',
                  'scripts/combine_tradis_plots',
                  'scripts/check_tradis_tags',
                  'scripts/remove_tradis_tags'])