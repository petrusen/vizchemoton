from setuptools import setup, find_packages

setup(
    name='vizchemoton',
    version='0.1.0',
    packages=find_packages(),
    install_requires=['networkx','bokeh>=2.3.2,<3.0.0',
                      'jsmol_bokeh_extension @ git+https://github.com/dgarayr/jsmol-bokeh-extension.git',
                      'matplotlib','numpy']
    python_requires='>=3.8',
    author=['Enric Petrus' 'Diego Garay-Ruiz'],
    author_email= ['enricpz@icloud.com'],
    description='A package for visualizing Chemoton data.',
    long_description=open('README.rst').read(),
    long_description_content_type='text/x-rst',
    url='https://github.com/petrusen/vizchemoton',
    include_package_data=True,
    keywords='chemistry visualization chemoton',
    project_urls={
        'Source': 'https://github.com/petrusen/vizchemoton',
        'Tracker': 'https://github.com/petrusen/vizchemoton/issues',
    },
)

