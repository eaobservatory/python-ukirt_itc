# Copyright (C) 2017 East Asian Observatory
# All Rights Reserved.
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>

from distutils.core import setup
import sys

sys.path.insert(0, 'lib')
from ukirt_itc.version import version

with open('README.rst') as f:
    long_description = f.read()

setup(
    name='ukirt_itc',
    version=version,
    description='UKIRT Integration Time Calculator',
    long_description=long_description,
    url='https://github.com/eaobservatory/python-ukirt_itc',
    package_dir={'': 'lib'},
    packages=['ukirt_itc'],
    package_data={'ukirt_itc': ['data/phot.json']})
