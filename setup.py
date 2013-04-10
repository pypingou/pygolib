#!/usr/bin/python
"""
Setup script
"""

from distutils.core import setup
from src import __version__

setup(
    name='pygolib',
    description='Tools to manipulate GO term and graphs.',
    author='Pierre-Yves Chibon',
    author_email='pingou@pingoured.fr',
    version=__version__,
    license='BSD',
    url='https://github.com/pypingou/pygolib/',
    package_dir={'pygolib': 'src'},
    packages=['pygolib'],
    scripts=["goutil"],
)
