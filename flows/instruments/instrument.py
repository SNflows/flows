from inspect import getmro

from astropy.time import Time
from astropy import units as u

from .. import FILTERS
from ..api import get_site

class MetaInstrument(type):

    def __new__(mcls, name, bases, attrs):

        all_attrs = dict()
        for base in bases:
            all_attrs.update(base.__dict__)
        all_attrs.update(attrs)

        if not 'siteid' in all_attrs:
            return super().__new__(mcls, name, bases, attrs)

        assert 'filters' in all_attrs, 'no filters variable'
        assert not set(all_attrs['filters'].values()) - FILTERS, 'unknown filter(s)'

        assert 'verify' in all_attrs, 'no verify classmethod'

        return super().__new__(mcls, name, bases, attrs)

class Instrument(metaclass=MetaInstrument):

    filter_keyword = 'filter'
    exptime_keyword = 'exptime'

    peakmax = None

    scale = None
    mirror = None
    angle = None

    def __init__(self):

        site = get_site(self.siteid)

        self.name = site['sitename']
        self.longitude = site['longitude']
        self.latitude = site['latitude']
        self.elevation = site['elevation']
        self.location = site['EarthLocation']

    def get_exptime(self, image):

        return float(image.hdr[self.exptime_keyword])

    def get_obstime(self, image):

        if 'mjd-obs' in image.hdr:
            obstime = Time(image.hdr['mjd-obs'], format='mjd', scale='utc', location=self.location)
        elif 'date-obs' in image.hdr:
            obstime = Time(image.hdr['date-obs'], format='isot', scale='utc', location=self.location)
        else:
            raise KeyError('mjd-obs nor date-obs in header')

        return obstime

    def get_filter(self, image):

        f = image.hdr[self.filter_keyword]

        assert f in self.filters, f'filter {f} not recognized'

        return self.filters[f]

    def get_mask(self, image):

        return None

    def __repr__(self):

        return self.name
