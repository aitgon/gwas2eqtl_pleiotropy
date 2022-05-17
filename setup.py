from setuptools import setup, find_packages

setup(
    name='eqtl2gwas_pleiotropy',
    version='0.1.0',
    url='https://tagc.univ-amu.fr/en/users/gonzalez-aitor',
    author='Aitor González',
    author_email='aitor.gonzalez@univ-amu.fr',
    description='Description of my package',
    packages=find_packages(),
    install_requires=['matplotlib', 'odfpy', 'pandas', 'requests', 'seaborn', 'sqlalchemy', 'termcolor'],
)
