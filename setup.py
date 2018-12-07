#!/usr/bin/env python
# -*- coding: utf-8 -*-
""" This script provides setup requirements to install hrvanalysis via pip"""

import setuptools

# Get long description in READ.md file
with open("README.md", "r") as fh:
    LONG_DESCRIPTION = fh.read()

setuptools.setup(
    name="hrv-analysis",
    version="1.0.3",
    author="Robin Champseix",
    license="GPLv3",
    author_email="robin.champseix@gmail.com",
    description="a package to calculate features from Rr Interval for HRV analyses",
    long_description=LONG_DESCRIPTION,
    long_description_content_type="text/markdown",
    include_package_data=True,
    url="https://github.com/robinchampseix/hrvanalysis",
    packages=setuptools.find_packages(),
    python_requires='>=3.5',
    install_requires=[
        "numpy>=1.15.1",
        "astropy>=3.0.4",
        "nolds>=0.4.1",
        "scipy>=1.1.0",
        "pandas>=0.23.4",
        "matplotlib>=2.2.2"
    ],
    classifiers=[
        "Programming Language :: Python :: 3",
        "Operating System :: OS Independent"
    ]
)
