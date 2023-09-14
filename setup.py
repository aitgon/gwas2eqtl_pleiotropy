from setuptools import setup, find_packages

setup(
    name='gwas2eqtl_pleiotropy',
    version='0.1.0',
    url='https://tagc.univ-amu.fr/en/users/gonzalez-aitor',
    author='Aitor González',
    author_email='aitor.gonzalez@univ-amu.fr',
    description='Description of my package',
    packages=find_packages(),
    install_requires=['crossmap', 'matplotlib', 'mysql-client', 'mysql-connector-python', 'odfpy', 'pandas',
                      'requests', 'seaborn>=0.12.2', 'snakemake', 'sqlalchemy', 'statannot', 'statsmodels',
                      'termcolor'],
)
