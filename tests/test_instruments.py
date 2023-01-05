"""
Test instruments module
"""
import logging

import pytest

from flows.instruments import INSTRUMENTS


def test_instruments():
    for instrument_name, instrument_class in INSTRUMENTS:
        instrument = instrument_class()
        site = instrument.get_site()
        logging.debug(f"{instrument_name}, site:{site['siteid']}, {instrument.siteid}")  # set log_cli_level=10 to show.
        if not instrument_name == "LCOGT":
            assert site['siteid'] == instrument.siteid
        else:
            assert site['siteid'] == instrument.image


if __name__ == '__main__':
    pytest.main([__file__])
