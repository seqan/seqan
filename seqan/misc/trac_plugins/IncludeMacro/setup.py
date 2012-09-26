#!/usr/bin/env python
# -*- coding: iso-8859-1 -*-
import os.path

from setuptools import setup

setup(
    name = 'TracIncludeMacro',
    version = '2.1fub3',
    packages = ['includemacro'],

    author = 'Manuel Holtgrewe',
    author_email = 'manuel.holtgrewe@fu-berlin.de',
    description = 'Include the contents of external URLs and other Trac objects in a wiki page.  Extended version that supports fragments.',
    long_description = open(os.path.join(os.path.dirname(__file__), 'README')).read(),
    license = 'BSD',
    keywords = 'trac 0.11 plugin wiki include macro',
    url = 'http://trac-hacks.org/wiki/IncludeMacro',
    classifiers = [
        'Framework :: Trac',
    ],
    
    install_requires = ['Trac'],

    entry_points = {
        'trac.plugins': [
            'includemacro.macros = includemacro.macros',
        ]
    }
)
