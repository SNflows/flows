import pickle, warnings

import numpy as np

from astropy.io import fits
from astropy.wcs import WCS, FITSFixedWarning

class Image:

    MASK_PARAMETERS = 'wcs',

    def __init__(self, data, hdr, wcs, exthdus=dict(), subtracted=None):

        self._data = data = np.asarray(data)
        self.x, self.y = np.meshgrid(*map(np.arange, data.shape[::-1]))

        self.hdr, self.wcs = hdr, wcs
        self.exthdus = exthdus

        self.set_subtracted(subtracted)

        self._lmasks = dict()
        self._mask, self._mask_hash = None, None
        self._data_mask = None

        self.instrument = None

    @classmethod
    def from_fits(cls, filename, subtracted=None):

        with fits.open(filename, mode='readonly') as hdul:

            data = hdul[0].data
            hdr = hdul[0].header
            exthdus = {hdu.name:hdu.copy() for hdu in hdul[1:]}

        with warnings.catch_warnings():

            warnings.simplefilter('ignore', category=FITSFixedWarning)
            wcs = WCS(hdr)

        if not subtracted is None:

            with fits.open(subtracted, mode='readonly') as hdul:
                subtracted = hdul[0].data

        else:

            subtracted = None

        return cls(data, hdr, wcs, exthdus, subtracted)

    def set_subtracted(self, subtracted):

        assert subtracted is None or np.shape(subtracted) == self._data.shape

        self._subtracted = np.asarray(subtracted) if not subtracted is None else None
        self._subtracted_mask = None

    def set_instrument(self, instrument):

        self.instrument = instrument
        self.filter = instrument.get_filter(self)
        self.obstime = instrument.get_obstime(self)
        self.exptime = instrument.get_exptime(self)

        if not (mask := instrument.get_mask(self)) is None:
            self.add_mask(mask)

    def add_mask(self, mask):

        i = max(self._lmasks.keys()) + 1 if len(self._lmasks) else 0
        self._lmasks[i] = mask

        return i

    def del_mask(self, i):

        assert i in self._lmasks, 'mask id does not exist'

        del self._lmasks[i]

    def _update_mask(self):

        mask_parameters = sorted([hash(mask) for mask in self._lmasks.values()])
        mask_parameters += [getattr(self, p) for p in self.MASK_PARAMETERS]
        mask_hash = hash(pickle.dumps(mask_parameters))

        print(self._mask_hash == mask_hash)
        if self._mask_hash == mask_hash:
            return False

        self._mask_hash = mask_hash

        masks = [~mask(self).astype(bool) for mask in self._lmasks.values()]
        self._mask = ~np.prod(masks, axis=0, dtype=bool) \
                if len(masks) else np.zeros_like(self._data, dtype=bool)

        return True

    @property
    def data(self):

        if self._update_mask() or self._data_mask is None:

            self._data_mask = self._mask.copy()
            self._data_mask |= ~np.isfinite(self._data)

        return np.ma.array(self._data, mask=self._data_mask)

    @property
    def subtracted(self):

        if self._subtracted is None:
            raise AttributeError('no subtracted image')

        if self._update_mask() or self._subtracted_mask is None:

            self._subtracted_mask = self._mask.copy()
            self._subtracted_mask |= ~np.isfinite(self._subtracted)

        return np.ma.array(self._subtracted, mask=self._subtracted_mask)

    def __getitem__(self, item):

        if item in ('filter', 'obstime', 'exptime') and \
                self.instrument is None:
            raise AttributeError('can not get %s without instrument' % item)

        return super().__getitem__(item)
