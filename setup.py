from setuptools import setup

requirements = [line.rstrip() for line in open("requirements.txt", "rt")]

with open('README.md') as f:
    readme = f.read()

setup(name='quatradis',
      version='0.1',
      description='Analyses transposon insertion sites',
      long_description=readme,
      url='http://quadram.ac.uk',
      author='Sarah Bastkowski',
      author_email='sarah.bastkowski@quadram.ac.uk',
      #license='MIT',   #TODO decide on license (original was GPL3)
      packages=['quatradis'],
      zip_safe=False,
      install_requires=requirements,
      scripts=[   'scripts/bacteria_tradis',
                  'scripts/tradis_plot',
                  'scripts/combine_tradis_plots',
                  'scripts/check_tradis_tags',
                  'scripts/remove_tradis_tags'])