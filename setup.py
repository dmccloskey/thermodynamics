#!/usr/bin/env python
# -*- coding: utf-8 -*-

from setuptools import setup

with open('README.rst') as readme_file:
    readme = readme_file.read()

with open('HISTORY.rst') as history_file:
    history = history_file.read()

requirements = [
    # TODO: put package requirements here
]

test_requirements = [
     # TODO: put package requirements here
]

setup(
    name='thermodynamics',
    version='0.1.0',
    description="thermodynamics module for cobrapy",
    long_description=readme + '\n\n' + history,
    author="Douglas McCloskey",
    author_email='domccl@biosustain.dtu.dk',
    url='https://github.com/Biosustain/thermodynamics',
    packages=[
        'thermodynamics',
    ],
    package_dir={'thermodynamics':
                 'thermodynamics'},
    include_package_data=True,
    install_requires=requirements,
    license="Apache Software License 2.0",
    zip_safe=False,
    keywords='thermodynamics',
    classifiers=[
        'Development Status :: 2 - Pre-Alpha',
        'Intended Audience :: Developers',
        'License :: OSI Approved :: Apache Software License',
        'Natural Language :: English',
        "Programming Language :: Python :: 2",
        'Programming Language :: Python :: 2.7',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.4',
        'Programming Language :: Python :: 3.5',
        'Programming Language :: Python :: 3.6',
    ],
    test_suite='tests',
    tests_require=test_requirements
)
