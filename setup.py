#!/usr/bin/env python
# -*- coding: utf-8 -*-

import setuptools
import hrvanalysis

# Get long description in READ.md file
with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="hrvanalysis",
    version=hrvanalysis.__version__,
    author="Robin Champseix",
    license="GPLv3",
    author_email="robin.champseix@gmail.com",
    description="a package to calculate features from Rr Interval for HRV analyses",
    long_description=long_description,
    long_description_content_type="text/markdown",
    include_package_data=True,
    url="https://github.com/robinchampseix/hrvanalysis",
    packages=setuptools.find_packages(),
    python_requires='>=3.6',
    setup_requires=[
        "numpy >= 1.15.1"
    ],
    install_requires=[
        "numpy >= 1.15.1",
        "astropy >= 3.0.4",
        "nolds >= 0.4.1",
        "scipy >= 1.1.0",
        "pandas >= 0.23.4",
        "matplotlib >= 2.2.2"
    ],
    classifiers=[
        "Programming Language :: Python :: 3",
        "Operating System :: OS Independent"
    ]
)
