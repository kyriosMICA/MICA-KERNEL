"""
MICA-Kernel v1.0
kyriosMICA · Benin, West Africa
TQIM-Davoh · Qudits-36
"""
from setuptools import setup, find_packages

setup(
    name            = 'mica-kernel',
    version         = '1.0.0',
    author          = 'Cyrille Egnon Davoh',
    author_email    = 'direction@kyriosmica.com',
    description     = 'TQIM-Davoh Analysis Engine — Qudits-36 Formalism',
    long_description= open('README.md').read(),
    long_description_content_type = 'text/markdown',
    url             = 'https://github.com/kyriosMICA/mica-kernel',
    packages        = find_packages(),
    python_requires = '>=3.9',
    install_requires= ['numpy>=1.21', 'scipy>=1.7'],
    classifiers     = [
        'Development Status :: 4 - Beta',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 3',
    ],
    entry_points    = {
        'console_scripts': ['mica=main:main'],
    },
)
