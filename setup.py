#!/usr/bin/env python
# -*- coding: utf-8 -*-


try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup


with open('README.rst') as readme_file:
    readme = readme_file.read()

with open('HISTORY.rst') as history_file:
    history = history_file.read().replace('.. :changelog:', '')

with open('requirements.txt') as requirements_file:
    requirements = requirements_file.read()

test_requirements = [
    'pytest', 'flake8'
]

setup(
    name='outrigger',
    version='0.2.8',
    description="Outrigger is a tool to de novo annotate splice sites "
                "and exons",
    long_description=readme + '\n\n' + history,
    author="Olga Botvinnik",
    author_email='olga.botvinnik@gmail.com',
    url='https://github.com/yeolab/outrigger',
    packages=[
        'outrigger', 'outrigger.index',
        'outrigger.psi', 'outrigger.io',
        'outrigger.validate'
    ],
    package_dir={'outrigger':
                 'outrigger'},
    include_package_data=True,
    install_requires=requirements,
    license="BSD",
    zip_safe=False,
    keywords='outrigger',
    classifiers=[
        'Development Status :: 2 - Pre-Alpha',
        'Intended Audience :: Developers',
        'License :: OSI Approved :: BSD License',
        'Natural Language :: English',
        "Programming Language :: Python :: 2",
        'Programming Language :: Python :: 2.7',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.4',
        'Programming Language :: Python :: 3.5',
    ],
    tests_require=test_requirements,
    entry_points={
        'console_scripts': [
            'outrigger = outrigger.commandline:main'
        ]
    },
)
