import os, logging

from inspect import getmro
from traceback import format_exc

from importlib import import_module

from .instrument import Instrument

INSTRUMENTS = list()

for instrument_file in os.listdir(__file__.rsplit('/',1)[0]):

    if not '.' in instrument_file:
        continue

    instrument_filename, file_ext = instrument_file.rsplit('.', 1)

    if instrument_filename[0] == '_' or file_ext != 'py':
        continue

    instrument = import_module('.' + instrument_filename, 'flows.instruments')

    for attribute in dir(instrument):

        if attribute[0] == '_':
            continue

        attribute = getattr(instrument, attribute)

        if not hasattr(attribute, '__bases__'):
            continue

        all_bases = [base for _base in attribute.__bases__ for base in getmro(_base)]
        if not Instrument in all_bases or not hasattr(attribute, 'siteid'):
            continue

        INSTRUMENTS.append(attribute)

def get_instrument(image):

    container = list()

    for instrument in INSTRUMENTS:
        try:
            instrument.verify(image)
        except Exception as e:
            if not str(e):
                e = format_exc().strip().split('\n')[-2].strip()
            logging.debug(f'{instrument} : {e}')
        else:
            container.append(instrument)

    assert container, 'Instrument was not identified'

    if len(container) > 1:
        msg = 'Data matched multiple instruments; '
        msg += ', '.join(map(str, container))
        logging.error(msg)

    return container[0]
