#!/usr/bin/env python3
import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="sage-projective-normgraph",
    version="0.1",
    author="Tomas Bayer",
    author_email="tomas.bayer@gmail.com",
    description="SageMath package for projective norm graphs",
    long_description=long_description,
    long_description_content_type="text/markdown",
    packages=setuptools.find_packages(),
    classifiers=[
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Mathematics",
        "Programming Language :: Python :: 2.7",
        "Operating System :: OS Independent"
    ]
)
