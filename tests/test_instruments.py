"""
Test instruments module
"""
import pytest
from flows.instruments import INSTRUMENTS


def test_instruments():
    for instrument_name, instrument_class in INSTRUMENTS:
        instrument = instrument_class()
        site = instrument.get_site()
        assert site['siteid'] == instrument.siteid


if __name__ == '__main__':
    pytest.main([__file__])
