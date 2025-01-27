from setuptools import setup, find_packages

setup(
    name='vizchemoton',
    version='0.1.0',
    packages=find_packages(),
    author='Enric Petrus, Diego Garay-Ruiz',
    author_email= ['enricpz@icloud.com'],
    description='A package for visualizing Chemoton data.',
    long_description=open('README.md').read(),
    url='https://github.com/petrusen/vizchemoton',
    license="BSD (3-clause)",
    include_package_data=True,
    keywords='chemistry visualization chemoton',
    project_urls={
        'Source': 'https://github.com/petrusen/vizchemoton',
        'Tracker': 'https://github.com/petrusen/vizchemoton/issues',
    },
)

