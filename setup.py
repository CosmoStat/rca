from setuptools import setup

setup(
    name='rca',
    version='2.0',
    description='Resolved Component Analysis',
    author='Morgan A. Schmitz, Fred Ngolè',
    author_email='morgan.schmitz@cea.fr',
    url='https://github.com/CosmoStat/rca',
    packages=['rca'],
    install_requires=['numpy','scipy','modopt']
)
