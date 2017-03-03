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
#
# This Python ITC is based on the original Perl UKIRT ITC.

from __future__ import absolute_import, division, print_function, \
    unicode_literals

from codecs import latin_1_decode
from collections import OrderedDict, namedtuple
import json
from math import exp, log10, pi, sqrt
from pkgutil import get_data

from .error import UKIRTITCError
from .version import version


SkyInfo = namedtuple(
    'SkyInfo',
    ('dark', 'grey', 'bright'))

InstrumentInfo = namedtuple(
    'InstrumentInfo',
    ('name', 'nread', 'dark', 'pixsize', 'gain', 'zeropoint'))

AreaInfo = namedtuple(
    'AreaInfo',
    ('pixel', 'object_', 'sky', 'frac_aperture', 'n_pix'))


class UKIRTImagPhotITC(object):
    # Instruments.
    UFTI = 1
    UIST = 2
    WFCAM = 3

    # Sky conditions.
    SKY_DARK = 1
    SKY_GREY = 2
    SKY_BRIGHT = 3

    # Limiting factors.
    LIMIT_BACKGROUND = 1
    LIMIT_READOUT = 2

    # Define valid filters and the order in which to show them.
    FILTERS = 'Z Y J H K L M BrG H2S1'.split()

    # Data structures to be filled with information from the data file.
    _sky = None
    _extinction = None
    _info = None

    def __init__(self):
        """
        Constructor method.
        """

        # Ensure we have read the data file.
        if self._info is None:
            self._read_data()

    @classmethod
    def _read_data(cls):
        """
        Read imaging photometry data from the file.
        """

        instrument_names = [
            (cls.UFTI, 'UFTI'),
            (cls.UIST, 'UIST'),
            (cls.WFCAM, 'WFCAM'),
        ]

        data = json.loads(latin_1_decode(
            get_data('ukirt_itc', 'data/phot.json'))[0])

        # Process sky data.
        data_sky = data.get('sky')
        if data_sky is None:
            raise UKIRTITCError('Data file did not contain "sky" section')
        cls._sky = {}
        for (filter_, values) in data_sky.items():
            if filter_ not in cls.FILTERS:
                raise UKIRTITCError(
                    'Sky filter "{0}" not recognised'.format(filter_))
            cls._sky[filter_] = SkyInfo(*values)

        # Process extinction data.
        data_ext = data.get('extinction')
        if data_ext is None:
            raise UKIRTITCError(
                'Data file did not contain "extinction" section')
        cls._extinction = {}
        for (filter_, value) in data_ext.items():
            if filter_ not in cls.FILTERS:
                raise UKIRTITCError(
                    'Extinction filter "{0}" not recognised'.format(filter_))
            cls._extinction[filter_] = value

        # Process instrument data.
        data_instruments = data.get('instrument')
        if data_instruments is None:
            raise UKIRTITCError(
                'Data file did not contain "instrument" section')
        cls._info = OrderedDict()
        for (instrument, name) in instrument_names:
            data_instrument = data_instruments.get(name)
            if data_instrument is None:
                raise UKIRTITCError(
                    'Could not find instrument information '
                    'for "{0}"'.format(name))

            info_obj = InstrumentInfo(name=name, **data_instrument)

            for filter_ in info_obj.zeropoint:
                if filter_ not in cls.FILTERS:
                    raise UKIRTITCError(
                        'Instrument "{0}" filter "{1}" not recognised'.format(
                            name, filter_))

            cls._info[instrument] = info_obj

    def get_version(self):
        """
        Get the module version.
        """

        return version

    def get_available_filters(self):
        """
        Get an ordered dictionary of filters, each giving a set of
        instruments for which that filter is available.
        """

        result = OrderedDict()

        for filter_ in self.FILTERS:
            avail = set()

            for (instrument, info) in self._info.items():
                if filter_ in info.zeropoint:
                    avail.add(instrument)

            if avail:
                result[filter_] = avail

        return result

    def calculate_time(
            self, mag, snr,
            instrument, filter_, aperture, sky, seeing, airmass,
            is_extended,
            with_extra_output=False):
        """
        Compute integration time for given magnitude and SNR.
        """

        instrument = self._info.get(instrument)
        if instrument is None:
            raise UKIRTITCError('Instrument not recognised.')

        zeropoint = instrument.zeropoint.get(filter_)
        if zeropoint is None:
            raise UKIRTITCError(
                'Filter "{0}" is not available for instrument "{1}".'.format(
                    filter_, instrument.name))

        area = self._determine_area(instrument, aperture, seeing, is_extended)

        if is_extended:
            mag = self._mag_per_pixel(mag, area)

        mag_sky = self._get_sky_magnitude(filter_, sky)

        extinction = self._extinction.get(filter_, 0.0)

        fd_sky = 10 ** (-0.4 * (mag_sky - zeropoint)) * instrument.gain  # Flux of the sky per e-/arcsec**2
        zp_airmass = zeropoint - (airmass - 1) * extinction
        flux_object = area.frac_aperture * 10 ** (-0.4 * (mag - zp_airmass)) * instrument.gain  # Flux of the source (e-/sec)

        a = (flux_object / snr) ** 2
        b = flux_object + fd_sky * area.object_ * (1 + area.object_ / area.sky) + instrument.dark * area.object_ / area.pixel * (1 + area.object_ / area.sky)
        c = instrument.nread ** 2 * area.object_ / area.pixel * (1 + area.object_ / area.sky)

        int_time = (b + sqrt(b ** 2 + 4 * a * c)) / (2 * a)

        if not with_extra_output:
            return int_time

        return (int_time, self._calculate_extra(
            instrument, area, int_time, flux_object, fd_sky))

    def calculate_magnitude(
            self, int_time, snr,
            instrument, filter_, aperture, sky, seeing, airmass,
            is_extended,
            with_extra_output=False):
        """
        Compute magnitude for given integraton time and SNR.
        """

        instrument = self._info.get(instrument)
        if instrument is None:
            raise UKIRTITCError('Instrument not recognised.')

        zeropoint = instrument.zeropoint.get(filter_)
        if zeropoint is None:
            raise UKIRTITCError(
                'Filter "{0}" is not available for instrument "{1}".'.format(
                    filter_, instrument.name))

        area = self._determine_area(instrument, aperture, seeing, is_extended)

        mag_sky = self._get_sky_magnitude(filter_, sky)

        extinction = self._extinction.get(filter_, 0.0)

        fd_sky = 10 ** (-0.4 * (mag_sky - zeropoint)) * instrument.gain  # Flux of the sky per e-/arcsec**2
        zp_airmass = zeropoint - (airmass-1) * extinction
        a = ((area.frac_aperture * int_time) ** 2) / (snr ** 2)
        b = -(int_time * area.frac_aperture)
        c = -(area.object_ / area.pixel) * (1 + area.object_ / area.sky) * (fd_sky * area.pixel * int_time + instrument.dark * int_time + instrument.nread ** 2)
        flux_object = (-1 * b + sqrt(b ** 2 - 4 * a * c)) / (2 * a)
        mag = -2.5 * log10(flux_object / instrument.gain) + zp_airmass
        flux_object = 10 ** (-0.4 * (mag - zp_airmass)) * instrument.gain * area.frac_aperture  # Flux of the source (e-/sec)

        if is_extended:
            mag = self._mag_per_sq_arcsec(mag, area)

        if not with_extra_output:
            return mag

        return (mag, self._calculate_extra(
            instrument, area, int_time, flux_object, fd_sky))

    def calculate_snr(
            self, int_time, mag,
            instrument, filter_, aperture, sky, seeing, airmass,
            is_extended,
            with_extra_output=False):
        """
        Compute SNR for given integration time and magnitude.
        """

        instrument = self._info.get(instrument)
        if instrument is None:
            raise UKIRTITCError('Instrument not recognised.')

        zeropoint = instrument.zeropoint.get(filter_)
        if zeropoint is None:
            raise UKIRTITCError(
                'Filter "{0}" is not available for instrument "{1}".'.format(
                    filter_, instrument.name))

        area = self._determine_area(instrument, aperture, seeing, is_extended)

        if is_extended:
            mag = self._mag_per_pixel(mag, area)

        mag_sky = self._get_sky_magnitude(filter_, sky)

        extinction = self._extinction.get(filter_, 0.0)

        fd_sky = 10 ** (-0.4 * (mag_sky - zeropoint)) * instrument.gain  # Flux of the sky per e-/arcsec**2
        zp_airmass = zeropoint - (airmass - 1) * extinction
        flux_object = area.frac_aperture * 10 ** (-0.4 * (mag - zp_airmass)) * instrument.gain  # Flux of the source (e-/sec)

        (nstar, nsky, ndark, noise) = self._calculate_noise(
            instrument, area, int_time, flux_object, fd_sky)

        snr = nstar / noise

        if not with_extra_output:
            return snr

        return (snr, self._calculate_extra(
            instrument, area, int_time, flux_object, fd_sky))

    def _determine_area(self, instrument, aperture, seeing, is_extended):
        """
        Compute instrument-specific area information.
        """

        area_pixel = instrument.pixsize ** 2

        if is_extended:
            # Aperture is 1 pixel and the sky aperture is the same.
            # There is no aperture loss.
            return AreaInfo(
                pixel=area_pixel, object_=area_pixel, sky=area_pixel,
                frac_aperture=1, n_pix=1)

        else:
            area_sky = 100
            area_object = pi * (aperture / 2) ** 2
            n_pix = int(area_object / area_pixel)

            # Fraction of flux inside aperture
            frac_aperture = 1 - exp(
                -(aperture / 2)**2 / (2 * (seeing / 2.35)**2))

            return AreaInfo(
                pixel=area_pixel, object_=area_object, sky=area_sky,
                frac_aperture=frac_aperture, n_pix=n_pix)

    def _calculate_noise(
            self, instrument, area, int_time, flux_object, fd_sky):
        """
        Compute noise information.
        """

        nstar = flux_object * int_time
        nsky = fd_sky * area.pixel * int_time
        ndark = instrument.dark * int_time
        noise = sqrt(nstar + (area.object_ / area.pixel) * (1 + area.object_ / area.sky) * (nsky + ndark + instrument.nread ** 2))

        return (nstar, nsky, ndark, noise)

    def _calculate_extra(
            self, instrument, area, int_time, flux_object, fd_sky):
        """
        Compute extra information to show with the results.
        """

        (nstar, nsky, ndark, noise) = self._calculate_noise(
            instrument, area, int_time, flux_object, fd_sky)

        return {
            'pixsize': instrument.pixsize,
            'n_pix': area.n_pix,
            'frac_aperture': area.frac_aperture,
            'e_obj': nstar,
            'e_sky': nsky,
            'noi_obj': sqrt(nstar),
            'noi_sky': sqrt((area.object_ / area.pixel * (1 + area.object_ / area.sky) * nsky)),
            'noi_ccd': sqrt((area.object_ / area.pixel * (1 + area.object_ / area.sky) * (ndark + instrument.nread ** 2))),
            'limit': (self.LIMIT_BACKGROUND if (nsky > (3 * ndark + instrument.nread ** 2)) else self.LIMIT_READOUT),
        }

    def _mag_per_pixel(self, mag, area):
        """
        Convert from mag/sq arcsec to mag/pixel.
        """

        return mag - 2.5 * log10(area.pixel)

    def _mag_per_sq_arcsec(self, mag, area):
        """
        Convert from mag/pixel to mag/sq arcsecond.
        """

        return mag + 2.5 * log10(area.pixel)

    def _get_sky_magnitude(self, filter_, sky):
        """
        Determine sky magnitude.
        """

        if filter_ not in self.FILTERS:
            raise UKIRTITCError('Filter {0} not recognised'.format(filter_))

        sky_info = self._sky.get(filter_)
        if sky_info is None:
            raise UKIRTITCError(
                'No sky information available for filter {0}'.format(filter_))

        if sky == self.SKY_DARK:
            return sky_info.dark

        if sky == self.SKY_GREY:
            return sky_info.grey

        if sky == self.SKY_BRIGHT:
            return sky_info.bright

        raise UKIRTITCError('Sky brightness not recognised')
