#!/usr/bin/python

import os
import setuptools

def readme():
    with open('README.md') as f:
        return f.read()

def get_requirements_filename():
        return "REQUIREMENTS.txt"

install_requires = [
    line.rstrip() for line in open(os.path.join(os.path.dirname(__file__), get_requirements_filename()))
]


setuptools.setup(
    name='antenna',
    version='0.0.4',
    description='A software package for identifying sgRNA reads in next-generation sequencing data',
    long_description=readme(),
    classifiers=[
      'Development Status :: 3 - Alpha',
      'Intended Audience :: Science/Research'
      'License :: OSI Approved :: BSD License',
      'Programming Language :: Python :: 3.7',
      'Topic :: Scientific/Engineering',
    ],
    url='http://github.com/broadinstitute/antenna',
    author='Nick Barkas',
    license='BSD (3-Clause)',
    packages=['antenna'],
    install_requires=install_requires,
    entry_points={
        'console_scripts': [
            'antenna_count_reads=antenna.antenna_count_reads:main',
            'antenna_tag_reads=antenna.antenna_tag_reads:main',     
        ],
    },
    include_package_data=True,
    zip_safe=False
)
