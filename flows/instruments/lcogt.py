from astropy.time import Time
from astropy.coordinates import EarthLocation

from .instrument import Instrument

class LCOGT(Instrument):

    filters = {
        'B' : 'B',
    }

    def get_obstime(self, image):

        lat, lon, height = image.hdr['latitude'], image.hdr['longitud'], image.hdr['height']
        location = EarthLocation.from_geodetic(lat, lon, height)

        return Time(image.hdr['mjd-obs'], format='mjd', scale='utc', location=location)

    def get_mask(self, image):

        return np.asarray(image.exthdus['BPM'].data, dtype=bool)

    @classmethod
    def verify(cls, image):

        assert image.hdr['origin'] == 'LCOGT'

class LCOGT_SAAO(LCOGT):

    siteid = 3

    @classmethod
    def verify(cls, image):

        LCOGT.verify(image)
        assert image.hdr['site'] == 'LCOGT node at SAAO'

class LCOGT_SSO(LCOGT):

    siteid = 6

    @classmethod
    def verify(cls, image):

        LCOGT.verify(image)
        assert image.hdr['site'] == 'LCOGT node at Siding Spring Observatory'
