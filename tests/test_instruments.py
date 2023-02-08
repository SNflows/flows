"""
Test instruments module
"""
import inspect
import logging

import pytest

import flows.instruments


def test_instruments():
    for instrument_name, instrument_class in flows.instruments.INSTRUMENTS:
        instrument = instrument_class()
        site = instrument.get_site()
        logging.debug(f"{instrument_name}, site:{site['siteid']}, {instrument.siteid}")  # set log_cli_level=10 to show.
        if not instrument_name == "LCOGT":
            assert site['siteid'] == instrument.siteid
        else:
            assert site['siteid'] == instrument.image


def test_all_instruments_present():
    ins = inspect.getmembers(flows.instruments.instruments,
                                 lambda member: inspect.isclass(member) and member.__module__ == flows.instruments.instruments.__name__)
    assert len(flows.instruments.INSTRUMENTS) == len(ins)

if __name__ == '__main__':
    pytest.main([__file__])
