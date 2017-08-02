#!/usr/bin/env python

try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup

setup(
    name='junction_counter',
    version='0.0.1',
    url='',
    license='',
    author='brianyee',
    author_email='',
    description='Count reads supporting intron inclusion or exclusion',
    packages=['junction_counter'],
    package_dir={
        'junction_counter': 'junction_counter',
    },
    entry_points = {
        'console_scripts': [
            'count-junctions = junction_counter.count_junctions:main',
        ]
    }
)