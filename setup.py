#!/usr/bin/env python
# -*- encoding: utf-8 -*-

from glob import glob
from os.path import basename, splitext

from setuptools import find_packages, setup

setup(
    name="snappynt",
    version="0.2.0",
    install_requires=["snappy>=3.0.1"],
    packages=find_packages("src"),
    package_dir={"": "src"},
    py_modules=[splitext(basename(path))[0] for path in glob("src/*.py")],
    package_data={"snappynt": ["data/*.json"]},
)
