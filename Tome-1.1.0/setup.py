#!/usr/bin/env python
# -*- coding: utf-8 -*-

from setuptools import setup, find_packages

long_description = open('README.md').read()
version = '1.1.0'

setup(
	name='tome',    # This is the name of your PyPI-package.
	description='Temperature optima for microorgianisms and enzymes',#package description
    long_description=long_description,
    version=version,                          # MAJOR.MINOR.PATCH
	author='Gang Li',
	author_email='gangl@chalmers.se',
	url='https://github.com/EngqvistLab/Tome',
    packages=find_packages(exclude=['test*']), #find folders containing scripts, exclude irrelevant ones
	# package_dir={'':'tome'},
    install_requires=[
        'pandas>=1.3.0',
        'Biopython>=1.78',
        'numpy>=1.21.0',
        'scikit-learn>=1.7.0',
        'scipy>=1.7.0',
        'joblib>=1.0.0',
        'requests>=2.25.0'
    ],
    include_package_data=True,
    package_data={'tome':['data/train.csv',
						  'model/OGT_svr.f',
						  'model/OGT_svr.pkl',
						  'external_data/enzyme_ogt_topt.tsv',
						  'external_data/brenda_sequences_20180109.fasta']},
	license='GPLv3+',
	classifiers=[
	# How mature is this project? Common values are
	#   3 - Alpha
	#   4 - Beta
	#   5 - Production/Stable
	'Development Status :: 3 - Alpha',

	# Indicate who your project is intended for
	'Intended Audience :: Science/Research',
	'Topic :: Scientific/Engineering',

	# Pick your license as you wish (should match "license" above)
	'License :: OSI Approved :: GNU General Public License v3 or later (GPLv3+)',

	# Specify the Python versions you support here
	'Programming Language :: Python :: 3',
	'Programming Language :: Python :: 3.8',
	'Programming Language :: Python :: 3.9',
	'Programming Language :: Python :: 3.10',
	'Programming Language :: Python :: 3.11',
	'Programming Language :: Python :: 3.12'],
    python_requires='>=3.8', #python version
    keywords='tome',

	entry_points={
        'console_scripts': ['tome = tome.tome:main']
    }
)
