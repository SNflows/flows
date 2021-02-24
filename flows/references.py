from collections import OrderedDict

import numpy as np

from astropy.table import Table, TableColumns, Column

class ReferenceColumns(TableColumns):

    def _set_image(self, image):

        self.image = image

    def keys(self):

        keys = list(super().keys())

        if not hasattr(self, 'image'):
            return keys

        keys += ['x'] if not 'x' in keys else []
        keys += ['y'] if not 'y' in keys else []

        return keys

    def values(self):

        return [self[key] for key in self.keys()]

    def __len__(self):

        return len(self.keys())

    def __iter__(self):

        return iter(key for key in self.keys())

    def __getitem__(self, item):

        if item in ('x', 'y') and \
                {'ra', 'decl'} <= set(self.keys()) and \
                hasattr(self, 'image'):

            rd = list(zip(self['ra'], self['decl']))
            x, y = zip(*self.image.wcs.all_world2pix(rd, 0))

            return Column({'x': x, 'y': y}[item], item)

        return super().__getitem__(item)

class References(Table):

    def __init__(self, *args, **kwargs):

        super().__init__(*args, **kwargs)

        self.columns = ReferenceColumns(self.columns)

    def set_image(self, image):

        self.columns._set_image(image)

    def __getitem__(self, item):

        if item == 'xy':
            return list(zip(self['x'], self['y']))

        return super().__getitem__(item)

    def copy_with_image(self):

        table = self.copy()
        del table['x'], table['y']
        table.set_image(self.columns.image)

        return table
