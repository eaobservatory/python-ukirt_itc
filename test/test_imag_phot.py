# Copyright (C) 2017 East Asian Observatory
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

from __future__ import absolute_import, division, print_function, \
    unicode_literals

from collections import OrderedDict
from sys import version_info
from unittest import TestCase

from ukirt_itc import UKIRTITCError, UKIRTImagPhotITC
from ukirt_itc.imag_phot import SkyInfo, InstrumentInfo

if version_info[0] < 3:
    string_type = unicode
else:
    string_type = str

BG = UKIRTImagPhotITC.LIMIT_BACKGROUND
RN = UKIRTImagPhotITC.LIMIT_READOUT


class ImagePhotTestCase(TestCase):
    def test_1_basic(self):
        """
        Basic test method.

        Note: must be run first (having method name which comes first
        by string ordering) as it tests filling of data structures.
        """

        # Error class should be an exception.
        self.assertIsInstance(UKIRTITCError(''), Exception)

        # Should start with empty data structures.
        self.assertIsNone(UKIRTImagPhotITC._sky)
        self.assertIsNone(UKIRTImagPhotITC._extinction)
        self.assertIsNone(UKIRTImagPhotITC._info)

        # After instantiating an ITC, data structures should have been filled.
        itc = UKIRTImagPhotITC()

        self.assertIsInstance(UKIRTImagPhotITC._sky, dict)
        for (key, value) in UKIRTImagPhotITC._sky.items():
            self.assertIsInstance(key, string_type)
            self.assertIn(key, UKIRTImagPhotITC.FILTERS)
            self.assertIsInstance(value, SkyInfo)

        self.assertIsInstance(UKIRTImagPhotITC._extinction, dict)
        for (key, value) in UKIRTImagPhotITC._extinction.items():
            self.assertIsInstance(key, string_type)
            self.assertIn(key, UKIRTImagPhotITC.FILTERS)
            self.assertIsInstance(value, float)

        self.assertIsInstance(UKIRTImagPhotITC._info, OrderedDict)
        for (key, value) in UKIRTImagPhotITC._info.items():
            self.assertIsInstance(key, int)
            self.assertIsInstance(value, InstrumentInfo)

    def test_2_wfcam(self):
        """
        Test WFCAM results against those from the original Perl ITC.
        """

        itc = UKIRTImagPhotITC()

        pixsize = 0.4

        param = {
            'instrument': UKIRTImagPhotITC.WFCAM,
            'aperture': 2.0,
            'sky': UKIRTImagPhotITC.SKY_GREY,
            'seeing': 0.9,
            'airmass': 1.2,
        }

        for (n, filter_, is_extended, int_time, mag, snr, n_pix, frac, e_obj, e_sky, noi_obj, noi_sky, noi_ccd, limit) in [
                #    Filter  Extended  time   mag  snr  npix  frac  e_obj    e_sky    obj     sky  ccd   lim
                (1,  'Z',    False,    21.0, 20.0, 5.0, 19.0, 0.97,  1194,    1137,  34.6,  151.8, 181.2, RN),
                (2,  'Y',    False,    29.0, 20.0, 5.0, 19.0, 0.97,  1560,    3105,  39.5,  250.8, 181.6, BG),
                (3,  'J',    False,    53.0, 20.0, 5.0, 19.0, 0.97,  3600,   23769,  60.0,  693.8, 183.0, BG),
                (4,  'H',    False,   259.0, 20.0, 5.0, 19.0, 0.97, 21230,  887319, 145.7, 4239.1, 194.1, BG),
                (5,  'K',    False,   767.0, 20.0, 5.0, 19.0, 0.97, 33721, 2242019, 183.6, 6738.3, 219.0, BG),
                (6,  'BrG',  False,  7842.0, 20.0, 5.0, 19.0, 0.97, 33520, 2208188, 183.1, 6687.3, 437.3, BG),
                (7,  'H2S1', False, 11891.0, 20.0, 5.0, 19.0, 0.97, 33581, 2212192, 183.3, 6693.3, 522.7, BG),

                (11, 'Z',    True,     51.0, 20.0, 5.0,  1.0, 1.0,    483,    2782,  22.0,   74.6,  57.5, BG),
                (12, 'Y',    True,     83.0, 20.0, 5.0,  1.0, 1.0,    739,    8895,  27.2,  133.4,  58.0, BG),
                (13, 'J',    True,    186.0, 20.0, 5.0,  1.0, 1.0,   2063,   82320,  45.4,  405.8,  59.8, BG),
                (14, 'H',    True,    936.0, 20.0, 5.0,  1.0, 1.0,  12664, 3198691, 112.5, 2529.3,  71.2, BG),
                (15, 'K',    True,   2767.0, 20.0, 5.0,  1.0, 1.0,  20125, 8086197, 141.9, 4021.5,  93.5, BG),
                (16, 'BrG',  True,  28284.0, 20.0, 5.0,  1.0, 1.0,  20005, 7964154, 141.4, 3991.0, 244.5, BG),
                (17, 'H2S1', True,  42888.0, 20.0, 5.0,  1.0, 1.0,  20041, 7978605, 141.6, 3994.6, 298.3, BG),

                ]:
            expect_extra = {
                'pixsize': pixsize,
                'n_pix': n_pix,
                'frac_aperture': frac,
                'e_obj': e_obj,
                'e_sky': e_sky,
                'noi_obj': noi_obj,
                'noi_sky': noi_sky,
                'noi_ccd': noi_ccd,
                'limit': limit,
            }

            # Test mode: calculate time.
            (result, extra) = itc.calculate_time(
                mag=mag, snr=snr, filter_=filter_,
                is_extended=is_extended, with_extra_output=True, **param)

            self.autoMsgAlmostEqual(result, int_time, n, 'time', 1.0)

            self._compare_extra(extra, expect_extra, n, 'time')

            # Test mode: calculate magnitude.
            (result, extra) = itc.calculate_magnitude(
                int_time=int_time, snr=snr, filter_=filter_,
                is_extended=is_extended, with_extra_output=True, **param)

            self.autoMsgAlmostEqual(result, mag, n, 'mag', 0.01)

            # TODO: self._compare_extra(extra, expect_extra, n, 'mag')

            # Test mode: calculate SNR.
            (result, extra) = itc.calculate_snr(
                int_time=int_time, mag=mag, filter_=filter_,
                is_extended=is_extended, with_extra_output=True, **param)

            self.autoMsgAlmostEqual(result, snr, n, 'snr', 0.1)

            # TODO: self._compare_extra(extra, expect_extra, n, 'mag')

    def _compare_extra(self, result, expect, n, t):
        self.autoMsgAlmostEqual(
            result['pixsize'], expect['pixsize'], n, t + ' pixsize', 0.001)

        self.autoMsgAlmostEqual(
            result['n_pix'], expect['n_pix'], n, t + ' n_pix', 0.001)

        self.autoMsgAlmostEqual(
            result['frac_aperture'], expect['frac_aperture'], n, t + ' frac. aperture', 0.01)

        self.autoMsgAlmostEqual(
            result['e_obj'], expect['e_obj'], n, t + ' e. object', 1)

        self.autoMsgAlmostEqual(
            result['e_sky'], expect['e_sky'], n, t + ' e. sky', 1)

        self.autoMsgAlmostEqual(
            result['noi_obj'], expect['noi_obj'], n, t + ' noise object', 0.1)

        self.autoMsgAlmostEqual(
            result['noi_sky'], expect['noi_sky'], n, t + ' noise sky', 0.1)

        self.autoMsgAlmostEqual(
            result['noi_ccd'], expect['noi_ccd'], n, t + ' noise CCD', 0.1)

        self.assertEqual(
            result['limit'], expect['limit'],
            'test #{0}: {1} limit got {2} expected {3}'.format(
                n, t, result['limit'], expect['limit']))

    def autoMsgAlmostEqual(self, value, expect, n, title, delta=None):
        kwargs = {}
        if delta is not None:
            kwargs['delta'] = delta

        self.assertAlmostEqual(
            value, expect,
            msg='test #{0}: {1}: got {2} expected {3}'.format(
                n, title, value, expect),
            **kwargs)
