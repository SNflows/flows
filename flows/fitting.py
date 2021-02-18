import numpy as np

from astropy.modeling.fitting import LevMarLSQFitter

class MaskableLevMarLSQFitter(LevMarLSQFitter):

    def __call__(self, model, x, y, z=None, *args, **kwargs):

        if hasattr(z, 'mask'):
            x, y, z = x[~z.mask], y[~z.mask], z[~z.mask]

        return super().__call__(model, x, y, z, *args, **kwargs)
