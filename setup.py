# setup.py
# Author: Xianjie Huang (xianjie5@connect.hku.hk)

from setuptools import setup, find_packages
from codecs import open
from os import path
import sys

here = path.abspath(path.dirname(__file__))

# load configures
exec(open("./xcltk/config.py").read())

# Get the long description from the relevant file
with open(path.join(here, "README.md"), encoding='utf-8') as f:
    long_description = f.read()

reqs = ['numpy>=1.9.0', 'pysam>=0.15.2', 'intervaltree', 'scipy', 'h5py',
        'pandas', 'anndata']

setup(
    name = "xcltk",

    # Versions should comply with PEP440.  For a discussion on single-sourcing
    # the version across setup.py and the project code, see
    # https://packaging.python.org/en/latest/single_source_version.html
    version = VERSION,

    description = "xcltk - Toolkit for XClone",
    long_description = long_description,
    long_description_content_type = "text/markdown",

    # The project's main homepage.
    url = "https://github.com/hxj5/xcltk",

    # Author details
    author = 'Xianjie Huang',
    author_email = 'xianjie5@connect.hku.hk',

    # Choose your license
    license='Apache-2.0',

    # What does your project relate to?
    keywords=['xclone', 'toolkit'],

    # You can just specify the packages manually here if your project is
    # simple. Or you can use find_packages().
    packages = find_packages(),

    #entry_points={
    #    'console_scripts': [
    #        'xcltk-baf = xcltk.baf.baf:main',
    #        'xcltk-rdr = xcltk.rdr.rdr:main',
    #        'xcltk-reg = xcltk.region.region:main'
    #    ],
    #}, 

    entry_points={
        'console_scripts': [
            'xcltk-plp-umi = xcltk.xcltk:main'
        ],
    }, 

    # List run-time dependencies here.  These will be installed by pip when
    # your project is installed. For an analysis of "install_requires" vs pip's
    # requirements files see:
    # https://packaging.python.org/en/latest/requirements.html
    
    install_requires = reqs,

    py_modules = ['xcltk']

    # buid the distribution: python setup.py sdist
    # upload to pypi: twine upload dist/...

)
