from astropy.time import Time
from astropy import units as u

from .instrument import Instrument

class Liverpool(Instrument):

    siteid = 8

    filters = {
        'Bessel-B' : 'B',
        'Bessell-B' : 'B',
        'V' : 'V' # XXX
    }

    filter_keyword = 'filter'#1' XXX

    def get_obstime(self, image):

        obstime = super().get_obstime(image)
        obstime += self.get_exptime(image) / 2 * u.second # Make time centre of exposure

        return obstime

    @classmethod
    def verify(cls, image):

        return
        assert image.hdr['telescop'] == 'Liverpool Telescope'
